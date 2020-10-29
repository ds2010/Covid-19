# Import pyomo module
from pyomo.environ import ConcreteModel, Set, Var, Objective, minimize, Constraint, log
from pyomo.opt import SolverFactory, SolverManagerFactory
from pyomo.core.expr.numvalue import NumericValue
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt


class CNLSZ:
    """Convex Nonparametric Least Square with contextual variables (CNLSZ)"""

    def __init__(self, y, x, z, cet='addi', fun='prod', rts='vrs'):
        """
            y : Output variable
            x : Input variables
            z : contexutal variables
            cet  = "addi" : Additive composite error term
                 = "mult" : Multiplicative composite error term
            fun  = "prod" : Production frontier
                 = "cost" : Cost frontier
            rts  = "vrs"  : Variable returns to scale
                 = "crs"  : Constant returns to scale
        """

        # TODO(error/warning handling): Check the configuration of the model exist
        self.x = x.tolist()
        self.y = y.tolist()    
        self.z = z.tolist() 
        self.cet = cet
        self.fun = fun
        self.rts = rts

        if type(self.x[0]) != list:
            self.x = []
            for x_value in x.tolist():
                self.x.append([x_value])
        
        if type(self.z[0]) != list:
            self.z = []
            for z_value in z.tolist():
                self.z.append([z_value])        

        # Initialize the CNLS model
        self.__model__ = ConcreteModel()

        # Initialize the sets
        self.__model__.I = Set(initialize=range(len(self.y)))
        self.__model__.J = Set(initialize=range(len(self.x[0])))
        self.__model__.K = Set(initialize=range(len(self.z[0])))

        # Initialize the variables
        self.__model__.alpha = Var(self.__model__.I, doc='alpha')
        self.__model__.beta = Var(self.__model__.I,
                                  self.__model__.J,
                                  bounds=(0.0, None),
                                  doc='beta')
        self.__model__.epsilon = Var(self.__model__.I, doc='residual')
        self.__model__.lamda = Var(self.__model__.K, doc='z coefficient')       
        self.__model__.frontier = Var(self.__model__.I,
                                      bounds=(0.0, None),
                                      doc='estimated frontier')

        # Setup the objective function and constraints
        self.__model__.objective = Objective(rule=self.__objective_rule(),
                                             sense=minimize,
                                             doc='objective function')
        self.__model__.regression_rule = Constraint(self.__model__.I,
                                                    rule=self.__regression_rule(),
                                                    doc='regression equation')
        if self.cet == "mult":
            self.__model__.log_rule = Constraint(self.__model__.I,
                                                 rule=self.__log_rule(),
                                                 doc='log-transformed regression equation')
        self.__model__.afriat_rule = Constraint(self.__model__.I,
                                                self.__model__.I,
                                                rule=self.__afriat_rule(),
                                                doc='afriat inequality')

        # Optimize model
        self.optimization_status = 0
        self.problem_status = 0

    def optimize(self, remote=True):
        """Optimize the function by requested method"""
        # TODO(error/warning handling): Check problem status after optimization
        if not remote:
            if self.cet == "addi":
                solver = SolverFactory("mosek")
                self.problem_status = solver.solve(self.__model__, tee=True)
                self.optimization_status = 1

            elif self.cet == "mult":
                solver = SolverFactory('gams')
                self.problem_status = solver.solve(self.__model__, solver='minos', tee=True)
                self.optimization_status = 1

        else:
            if self.cet == "addi":
                opt = "mosek"

            elif self.cet == "mult":
                opt = "knitro"

            solver = SolverManagerFactory('neos')
            self.problem_status = solver.solve(self.__model__,
                                               tee=True,
                                               opt=opt)
            self.optimization_status = 1

    def __objective_rule(self):
        """Return the proper objective function"""

        def objective_rule(model):
            return sum(model.epsilon[i] ** 2 for i in model.I)

        return objective_rule

    def __regression_rule(self):
        """Return the proper regression constraint"""
        if self.cet == "addi":
            if self.rts == "vrs":

                def regression_rule(model, i):
                    return self.y[i] == model.alpha[i] \
                        + sum(model.beta[i, j] * self.x[i][j] for j in model.J) \
                        + sum(model.lamda[k] * self.z[i][k] for k in model.K) + model.epsilon[i]

                return regression_rule
            elif self.rts == "crs":
                # TODO(warning handling): replace with model requested not exist
                return False

        elif self.cet == "mult":

            def regression_rule(model, i):
                return log(self.y[i]) == log(model.frontier[i] + 1) \
                        + sum(model.lamda[k] * self.z[i][k] for k in model.K) + model.epsilon[i]

            return regression_rule

        # TODO(error handling): replace with undefined model attribute
        return False

    def __log_rule(self):
        """Return the proper log constraint"""
        if self.cet == "mult":
            if self.rts == "vrs":

                def log_rule(model, i):
                    return model.frontier[i] == model.alpha[i] + sum(
                        model.beta[i, j] * self.x[i][j] for j in model.J) - 1

                return log_rule
            elif self.rts == "crs":

                def log_rule(model, i):
                    return model.frontier[i] == sum(
                        model.beta[i, j] * self.x[i][j] for j in model.J) - 1

                return log_rule

        # TODO(error handling): replace with undefined model attribute
        return False

    def __afriat_rule(self):
        """Return the proper afriat inequality constraint"""
        if self.fun == "prod":
            __operator = NumericValue.__le__
        elif self.fun == "cost":
            __operator = NumericValue.__ge__

        if self.cet == "addi":
            if self.rts == "vrs":

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                             for j in model.J),
                        model.alpha[h] + sum(model.beta[h, j] * self.x[i][j]
                                             for j in model.J))

                return afriat_rule
            elif self.rts == "crs":
                # TODO(warning handling): replace with model requested not exist
                return False
        elif self.cet == "mult":
            if self.rts == "vrs":

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        model.alpha[i] + sum(model.beta[i, j] * self.x[i][j]
                                             for j in model.J),
                        model.alpha[h] + sum(model.beta[h, j] * self.x[i][j]
                                             for j in model.J))

                return afriat_rule
            elif self.rts == "crs":

                def afriat_rule(model, i, h):
                    if i == h:
                        return Constraint.Skip
                    return __operator(
                        sum(model.beta[i, j] * self.x[i][j] for j in model.J),
                        sum(model.beta[h, j] * self.x[i][j] for j in model.J))

                return afriat_rule

        # TODO(error handling): replace with undefined model attribute
        return False

    def display_status(self):
        """Display the status of problem"""
        if self.optimization_status == 0:
            self.optimize()
        print(self.display_status)

    def display_alpha(self):
        """Display alpha value"""
        if self.optimization_status == 0:
            self.optimize()
        self.__model__.alpha.display()

    def display_beta(self):
        """Display beta value"""
        if self.optimization_status == 0:
            self.optimize()
        self.__model__.beta.display()

    def display_residual(self):
        """Dispaly residual value"""
        if self.optimization_status == 0:
            self.optimize()
        self.__model__.epsilon.display()

    def get_status(self):
        """Return status"""
        return self.optimization_status

    def get_alpha(self):
        """Return alpha value by array"""
        if self.optimization_status == 0:
            self.optimize()
        alpha = list(self.__model__.alpha[:].value)
        return np.asarray(alpha)

    def get_beta(self):
        """Return beta value by array"""
        if self.optimization_status == 0:
            self.optimize()
        beta = np.asarray([i + tuple([j]) for i, j in zip(list(self.__model__.beta),
                                                          list(self.__model__.beta[:, :].value))])
        beta = pd.DataFrame(beta, columns=['Name', 'Key', 'Value'])
        beta = beta.pivot(index='Name', columns='Key', values='Value')
        return beta.to_numpy()

    def get_lamda(self):
        """Return beta value by array"""
        if self.optimization_status == 0:
            self.optimize()
        lamda = list(self.__model__.lamda[:].value)
        return np.asarray(lamda)
    
    def get_residual(self):
        """Return residual value by array"""
        if self.optimization_status == 0:
            self.optimize()
        residual = list(self.__model__.epsilon[:].value)
        return np.asarray(residual)

    def get_frontier(self):
        """Return estimated frontier value by array"""
        if self.optimization_status == 0:
            self.optimize()
        frontier = list(self.__model__.frontier[:].value)
        return np.asarray(frontier)

    def get_adjusted_residual(self):
        """Return the shifted residuals(epsilon) tern by CCNLS"""
        return self.get_residual() - np.amax(self.get_residual())

    def get_adjusted_alpha(self):
        """Return the shifted constatnt(alpha) term by CCNLS"""
        return self.get_alpha() + np.amax(self.get_residual())

    def plot2d(self, xselect, fig_name=None):
        """Plot with selected x"""
        x = np.array(self.x).T[xselect]
        y = np.array(self.y).T
        f = y - self.get_residual()
        data = (np.stack([x, y, f], axis=0)).T

        # sort
        data = data[np.argsort(data[:, 0])].T

        x, y, f = data[0], data[1], data[2]

        # create figure and axes objects
        fig, ax = plt.subplots()
        dp = ax.scatter(x, y, color="k", marker='x')
        fl = ax.plot(x, f, color="r", label="CNLS")

        # add legend
        legend = plt.legend([dp, fl[0]],
                            ['Data points', 'CNLS'],
                            loc='upper left',
                            ncol=1,
                            fontsize=10,
                            frameon=False)

        # add x, y label
        ax.set_xlabel("Input $x$%d" % (xselect))
        ax.set_ylabel("Output $y$")

        # Remove top and right axes
        ax.spines["right"].set_visible(False)
        ax.spines["top"].set_visible(False)
        if fig_name == None:
            plt.show()
        else:
            plt.savefig(fig_name)
            
