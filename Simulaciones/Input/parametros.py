def parametros_ODE():
    eq = input("¿Que ecuacion vas a modelar?")
    parametros = []
    initial_condition = float(input('Ingrese una condicion inicial:'))
    if eq == "logistica":
        k = float(input('Ingrese k:'))
        r = float(input('Ingrese r:'))
        parametros.append(k)
        parametros.append(r)
    return eq, parametros, initial_condition


def parametros_PDE(eq, fuente):
    if eq == 'pndls':
        parametros = [fuente, 1, 1, 1, 0.28, 0.1, 0.32] #[forzamiento, control, alpha, beta, gamma, mu, nu]
    elif eq == 'wave':
        parametros = [5] #[c]
    elif eq == 'fisher':
        parametros = [5, 0.1]
    elif eq == 'KdV':
        parametros = [1, 1]
    return parametros
