from dtk.utils.analyzers.BaseAnalyzer import BaseAnalyzer
from dtk.utils.analyzers.timeseries import TimeseriesAnalyzer
from relative_time import *
from simtools.AnalyzeManager.AnalyzeManager import AnalyzeManager
from simtools.SetupParser import SetupParser
from simtools.Utilities.Experiments import retrieve_experiment
import numpy as np
import pandas as pd


def new_plot_fn(df, ax):
    grouped = df.groupby(level=['group'], axis=1)
    m = grouped.mean()
    m.plot(ax=ax,legend=True)

def new_filter_fn(sim_metadata):
    # return (sim_metadata['Run_Number'] == 0)
    if isinstance(sim_metadata['IRS'], bool):
        all_intervene_on = (
        sim_metadata['IRS'] and sim_metadata['ITNs'] and sim_metadata['MDA'] and sim_metadata['MSAT'] and sim_metadata[
            'StepD'] and sim_metadata['Healthseek'])
    else:
        all_intervene_on = (sim_metadata['IRS'] == "true" and
                            sim_metadata['ITNs'] == "true" and
                            sim_metadata['MDA'] == "true" and
                            sim_metadata['MSAT'] == "true" and
                            sim_metadata['StepD'] == "true" and
                            sim_metadata['Healthseek'] == "true")
    return all_intervene_on and (sim_metadata['Run_Number'] == 0)

if __name__=="__main__":
    SetupParser.init('HPC')

    am = AnalyzeManager()

    am.add_experiment(retrieve_experiment("c4f143ce-d60c-e711-9400-f0921c16849c")) #Cait
    # am.add_experiment(retrieve_experiment("40221d78-b289-e711-9401-f0921c16849d"))  # old .exe with new climate
    # am.add_experiment(retrieve_experiment("5e854ab1-bb89-e711-9401-f0921c16849d")) # old .exe with old climate

    # am.add_experiment(retrieve_experiment("40221d78-b289-e711-9401-f0921c16849d"))  # old .exe with old climate and enable_demographics_other = 0
    # am.add_experiment(retrieve_experiment("5e854ab1-bb89-e711-9401-f0921c16849d")) # old .exe with old climate and enable_demographics_other = 1

    # am.add_experiment(retrieve_experiment("4a2a24b5-f589-e711-9401-f0921c16849d"))  # deconstruct_L1
    # am.add_experiment(retrieve_experiment("f2e1a7e7-308c-e711-9401-f0921c16849d"))  # deconstruct_L4
    # am.add_experiment(retrieve_experiment("001a9f44-758c-e711-9401-f0921c16849d")) # deconstruct_L5
    # am.add_experiment(retrieve_experiment("4188b9de-e28c-e711-9401-f0921c16849d"))  # L6
    # am.add_experiment(retrieve_experiment("d939e2d9-fc8c-e711-9401-f0921c16849d"))  # new L6_pop1000_fixedstepd

    # am.add_experiment(retrieve_experiment("bbbfdbdf-fc8c-e711-9401-f0921c16849d"))
    # am.add_experiment(retrieve_experiment("d939e2d9-fc8c-e711-9401-f0921c16849d")) #1000
    # # am.add_experiment(retrieve_experiment("d839e2d9-fc8c-e711-9401-f0921c16849d"))
    # am.add_experiment(retrieve_experiment("ed60b5d2-fc8c-e711-9401-f0921c16849d")) #4761

    am.add_experiment(retrieve_experiment("2007274b-3c8d-e711-9401-f0921c16849d"))  # L10 4761
    am.add_experiment(retrieve_experiment("fd452d45-3c8d-e711-9401-f0921c16849d"))  # L10 1000

    am.add_analyzer(TimeseriesAnalyzer(plot_function=new_plot_fn,
                                       filter_function=new_filter_fn,
                                       channels=['Rainfall', 'Adult Vectors','Daily Bites per Human','New Diagnostic Prevalence','Air Temperature']))
    am.analyze()