import src as zn

sim = zn.Sim()
new_pars = {'rand_seed': 1}
sim.update_pars(new_pars)
msim = zn.MultiSim(sim, n_runs = 10)

if __name__ == "__main__":
    msim.run(verbose = 2, parallel = True, n_cpus = 4)
