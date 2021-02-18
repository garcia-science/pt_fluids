def parametros():
    eq = input("¿Que ecuacion vas a modelar?")
    parametros = []
    initial_condition = float(input('Ingrese una condicion inicial:'))
    if eq == "logistica":
        k = float(input('Ingrese k:'))
        r = float(input('Ingrese r:'))
        parametros.append(k)
        parametros.append(r)
    return eq, parametros, initial_condition