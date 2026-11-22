# Working on the experiments

Keep paper presets explicit in the configuration and suite functions. When changing a coefficient, record the model, backward variant, time step, window and observation mask together. Use a new results directory and retain the saved settings and source hashes.

Run tests/verify_numerics after numerical changes and tests/verify_resume after changes to run reuse or data identity. Use the smoke suites when changing output or plotting. Reproduce an affected full experiment when its numerical behavior changes. Failed runs must remain labeled as failed and must not acquire completed-cycle endpoints.

Reference figures are published experiment snapshots. Replace them only together with the associated numerical data, figure manifest and validation record. Keep supplemental archive plots distinct from the manuscript reference figures.
