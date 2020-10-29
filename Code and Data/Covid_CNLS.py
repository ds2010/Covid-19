import pandas as pd
import numpy as np
import CNLSZ
import time
from pyomo.environ import Constraint


# read covid-19 data
df = pd.read_excel('covid.xlsx', index_col=None)

# output
y = df['Dead']

# inputs
x1 = df['Covid Bed Occupancy']
x1 = np.asmatrix(x1).T
x2 = df['Covid MV bed occupancy']
x2 = np.asmatrix(x2).T
x = np.concatenate((x1, x2), axis=1)

# contextual variables
z1 = df['Share of 65-64']
z1 = np.asmatrix(z1).T
z2 = df['Share of +85']
z2 = np.asmatrix(z2).T
z3 = df['Staff absence / max weekly bed occupancy']
z3 = np.asmatrix(z3).T
z4 = df['week']
z4 = np.asmatrix(z4).T

z = np.concatenate((z1, z2, z3, z4), axis=1)

# define and solve the CNLS model
start = time.time()
res = CNLSZ.CNLSZ(y, x, z, cet= "mult", fun= "prod", rts= "vrs")

def constraint_rule(model, i, j):
    upperbound = [1, 1]
    return model.beta[i, j] <= upperbound[j]

res.__model__.beta_constraint_rule = Constraint(res.__model__.I,
                                                res.__model__.J,
                                                rule=constraint_rule,
                                                doc='beta constraint')
res.optimize(remote=False)
print("Time used :", time.time() - start)

# convert results to DataFrame
alpha = pd.DataFrame(res.get_alpha())
alpha.columns = ['Intercept']
beta = pd.DataFrame(res.get_beta())
beta.columns = ['Covid Bed Occupancy', 'Covid MV bed occupancy']
lamda = pd.DataFrame(res.get_lamda())
lamda.columns = ['Contextual variables']
residual = pd.DataFrame(res.get_residual())
residual.columns = ['residuals']

# Create a Pandas Excel writer using XlsxWriter as the engine.
writer = pd.ExcelWriter('vrs_covid_results.xlsx', engine='xlsxwriter')

# Write each dataframe to a different worksheet.
alpha.to_excel(writer, sheet_name='alpha')
beta.to_excel(writer, sheet_name='beta')
lamda.to_excel(writer, sheet_name='lamda')
residual.to_excel(writer, sheet_name='residual')

# Close the Pandas Excel writer and output the Excel file.
writer.save()


