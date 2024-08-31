import pandas as pd
from datetime import timedelta
import numpy as np
import pystan
import pickle

def get_stan_model_code():
    stan_model_code = """
        data {{
            int<lower = 0> m; // number of time points (T in the paper notation)
            int<lower = 0> k; // number of latent factors
            int<lower = 0> n; // number of counties
            int<lower=0, upper=1> y_miss[n]; //indicator of missingness in the final col of Y
            int<lower=0> y[n,m]; //Y(0) matrix with counties in rows and times in columns
            int<lower=0> offs[n]; //offset (population size)
        }}

        transformed data {{
            int<lower = 1> nmiss; //number of missings

            nmiss=sum(y_miss);
        }}

        parameters {{
            real alpha; //Global intercept
            vector[n] d0; //County deviations from global intercept
            row_vector[m] c0; //time-specific deviations from global intercept
            matrix[n,k] FS; //Factor scores matrix (U in the grant notation)
            matrix[m, k] L; //factor loadings matrix (V in the grant notation)
            real<lower=0> phi; //neg binomial scale parameter
        }}

        transformed parameters {{
            matrix[n,m] Ups; //intermediate predictor
            row_vector<lower=0>[m] Mu[n]; //negative binomial mean
            
            // latent predictors
            Ups = FS * L';
            // mean structure
            //for(i in 1:n) Mu[i] = exp(alpha + d0[i] + c0 + Ups[i] + log(offs)[i]);
            for (i in 1:n) {{
                for (j in 1:m) {{
                    real log_mu_pred = alpha + d0[i] + c0[j] + Ups[i, j] + log(offs[i]);
                    if (log_mu_pred > 15) {{
                        Mu[i, j] = exp(15);
                    }} else {{
                        Mu[i, j] = exp(log_mu_pred);
                    }}
                }}
            }}
        }}

        model {{

            // Priors for parameters
            alpha ~ normal(0, 5); // Prior for global intercept
            d0 ~ normal(0, 5); // Prior for county deviations
            c0 ~ normal(0, 5); // Prior for time deviations
            to_vector(FS) ~ normal(0, 1); // Prior
            to_vector(L) ~ normal(0, 1); // Prior
            phi ~ lognormal(0, 1); // Log-normal prior for positive scale parameter


            for(i in 1:n){{
                if (1-y_miss[i]) y[i] ~ neg_binomial_2(Mu[i],phi); //Likelihood contribution for control counties
                else y[i,1:(m-5)] ~ neg_binomial_2(Mu[i,1:(m-5)],phi); //Likelihood contribution for treated counties
            }}
        }}

        generated quantities {{
            int<lower=0> Y_pred[nmiss, 5]; // Compute the predictions for treated units at treated times
            {{
                int idy = 0;
                for (i in 1:n) {{
                    if (y_miss[i]) {{
                        for (j in 1:5) {{
                            Y_pred[idy + 1, j] = neg_binomial_2_rng(Mu[i, m - 5 + j], phi); // Predict last 5 columns
                        }}
                        idy = idy + 1;
                    }}
                }}
            }}
        }}
    """
    return stan_model_code

def get_population(state, year, pop_df):
    state_pop = pop_df.loc[pop_df['ST_NM'] == state, pop_df['year'] == year].copy()
    return state_pop['population']

def get_offset(disease_matrix, id, pop_df, cyclone_df):
    cyclone_data = cyclone_df[cyclone_df['sid'] == id]
    year = pd.to_datetime(cyclone_data['start_time'].values[0]).year
    pop_list = []
    for state in disease_matrix.index:
        population = get_population(state, year, pop_df)
        pop_list.append(population)

    return np.array(pop_list)

def run_indiv_model(tc_id, stan_model, panel_matrix, disease_matrix, pop_df, cyclone_df, k = 4, adapt_delta = 0.99, max_treedepth = 15):
    print("running model for TC ID:", tc_id)
    
    y_miss = panel_matrix.iloc[:, -1] == 1
    y_miss_numeric = y_miss.astype(int)
    y = disease_matrix.astype(float)
    y.loc[y_miss, y.columns[-5:]] = np.nan
    y[np.isnan(y)] = 99999
    offset = get_offset(disease_matrix, tc_id, pop_df, cyclone_df)

    stan_data = {
        'k': k,
        'm': y.shape[1],
        'n': y.shape[0],
        'y_miss': y_miss_numeric,
        'y': y.astype(int),
        'offs': offset.astype(int)
        }
    
    try:
        fit = stan_model.sampling(data=stan_data, chains=2, iter=2000, warmup=1000, 
                                  control={'adapt_delta': adapt_delta, 'max_treedepth': max_treedepth})
        
        samples = fit.extract(permuted=True)['Y_pred']
        hmc_diagnostics = pystan.diagnostics.check_hmc_diagnostics(fit)
        summary_dict = fit.summary()
        rhat_values = summary_dict['summary'][:, -1]
        param_names = summary_dict['summary_rownames']
        rhat_dict = dict(zip(param_names, rhat_values))

        return samples, hmc_diagnostics, rhat_dict

    except RuntimeError as e:
        print(f"Error running model for TC ID {tc_id}: {e}")
        return None
    
def run_submodels_state(panel_matrices, disease_matrices, stan_model, pop_df, cyclone_df, pickle_file, k = 4, adapt_delta = 0.99, max_treedepth = 15):
    Y0_psamp_dict = {}
    fit_dict = {}

    for tc_id in panel_matrices.keys():
        panel_matrix = panel_matrices.get(tc_id)
        disease_matrix = disease_matrices.get(tc_id)
        result = run_indiv_model(tc_id, stan_model, panel_matrix, disease_matrix, pop_df, cyclone_df, k, adapt_delta, max_treedepth)
        if result is not None:
            samples, diagnostics, r_hat_dict = result
            fit_dict[tc_id] = {'diagnostics': diagnostics, 'r_hat': r_hat_dict}
            Y0_psamp_dict[tc_id] = samples

    with open(pickle_file, 'wb') as f:
        pickle.dump({'samples': Y0_psamp_dict, 'fit_data': fit_dict}, f)