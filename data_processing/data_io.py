import datetime
import os
import pandas as pd

def read_rga_outs(filepath):
    """Reads in all '.rga_out' files in filepath as RGA data output files and returns the data as a pandas DataFrame.
    
    Keyword arguments:
    filepath -- the path to the folder containing RGA text files
    """
    files = os.listdir(filepath)
    files = [value for value in files if value.endswith(".rga_out")]
    blocks = []

    for filename in files:

        rga_times = []
        atomic_masses = []
        partial_pressures = []
        with open(filename, 'r') as f:

            for idx, line in enumerate(f):

                if idx == 0:
                    rga_time = datetime.datetime.strptime(line.rstrip("\n"), "%B %d, %Y  %I:%M:%S %p")

                elif idx > 21:
                    line = line.replace(" ", "").rstrip(",\n")
                    atomic_mass, partial_pressure = line.split(",")
                    rga_times.append(rga_time)
                    atomic_masses.append(float(atomic_mass))
                    partial_pressures.append(float(partial_pressure))

        data_block = pd.DataFrame(data=[rga_times, atomic_masses, partial_pressures]).T
        data_block.columns = ['Datetime', 'Mass_amu', 'Partial_Pressure_torr']
        blocks.append(data_block)

    return pd.concat(blocks)
