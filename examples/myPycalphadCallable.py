from pycalphad import Database, equilibrium, variables as v
from pycalphad.core.utils import instantiate_models, filter_phases, unpack_components
from pycalphad.codegen.callables import build_phase_records

# Problem Specific setup for our 9-element space exploration. Make sure that the elements 
# are in the same order as the composition vector that will be passed to the equilibrium_callable.
dbf = Database("CrHfMoNbTaTiVWZr_9element_Feb2023.tdb")
T = 1000
elementalSpaceComponents = ['Ti', 'Zr', 'Hf', 'W', 'Nb', 'Ta', 'Mo']

# Setup the pycalphad models
phases = list(set(dbf.phases.keys()))
comps = [s.upper() for s in elementalSpaceComponents]+['VA']
phases_filtered = filter_phases(dbf, unpack_components(dbf, comps), phases)
models = instantiate_models(dbf, comps, phases_filtered)
phase_records = build_phase_records(dbf, comps, phases_filtered, {v.N, v.P, v.T}, models=models)
expected_conds=[v.T]+[v.X(el) for el in comps[:-2]]
default_conds={v.P: 101325, v.N: 1.0}

# A neat callable for the equilibrium calculation that we will pass to the parallel graph exploration
def equilibrium_callable(elP):
    # Round to 6 decimal places, but make sure that 0.0 is not rounded to 0.0. pyCalphad will not like 
    # if the components are exactly 0 or sum to exactly 1 in the equilibrium calculation for reasons
    # that are beyond the scope of this tutorial.
    elP_round = [round(v-0.000001, 6) if v>0.000001 else 0.0000001 for v in elP]
    conds = {**default_conds, **dict(zip(expected_conds, [T] + elP_round))}
    eq_res = equilibrium(
        dbf, comps, phases_filtered, 
        conds, model=models, phase_records=phase_records, calc_opts=dict(pdens=2000))
    # Process the result into what we want -> a list of phases present
    nPhases = eq_res.Phase.data.shape[-1]
    phaseList = list(eq_res.Phase.data.reshape([nPhases]))
    phasePresentList = [pn for pn in phaseList if pn!='']
    if len(phasePresentList)==0:
        # Decide on what to do if the equilibrium calculation failed to converge. By default, we will
        # just pass and let the graph exploration scheme decide what to do about it.
        # print(f"Point: {elP_round} failed to converge.")
        pass
    return phasePresentList

# # Some extra code for the future tutorial on Scheil solidification :)
# from scheil import simulate_scheil_solidification
# # Meta Settings
# scheil_start_temperature = 3000
# liquid_phase_name = 'LIQUID'
# def scheil_callable(elP):
#     elP_round = [round(v-0.000001, 6) if v>0.000001 else 0.0000001 for v in elP]
#     initial_composition = {**dict(zip([v.X(el) for el in comps[:-2]], elP_round))}
#     
#     sol_res = simulate_scheil_solidification(
#         dbf, comps, phases_filtered, 
#         initial_composition, scheil_start_temperature, step_temperature=1.0)
# 
#     phaseFractions = {}
#     for phase, ammounts in sol_res.cum_phase_amounts.items():
#         finalAmmount = round(ammounts[-1], 6)
#         if finalAmmount>0:
#             phaseFractions.update({phase: finalAmmount})
# 
#     return phaseFractions.keys()

if __name__ == "main":
    pass