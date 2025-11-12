import zoonosim as zn

new_pars = dict(
    start_day = '2022-01-01',
    end_day = '2025-12-31'
)

sim = zn.Sim(datafile="zoonosim/data/H5N1_cases_in_QC_poultry.csv", pars=new_pars)
sim.run()
sim.compute_fit()
print(sim.fit.mismatch)
# fit.plot()
