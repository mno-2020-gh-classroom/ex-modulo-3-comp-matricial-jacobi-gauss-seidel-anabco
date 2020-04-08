########################## Información general ##########################
# El método Jacobi permite resolver sistemas de ecuaciones del tipo Ax=b.
# Esto se logra mediante un proceso iterativo.
# En esta primer versión, el número de iteraciones es fijo

########################## Declaración de funciones ##########################
# Función que obtiene cada componente del vector
funcObtenerComponente <- function(i, n, mtrx_A, vct_X, vct_B){

  # Variable en la cual acumularemos el resutlado de la sumatoria
  nbr_Sumatoria = 0
  nbr_Final = 0

  # Variables con los términos agrupados de los elementos de la fórmula
  nbr_Termino1 = 0
  nbr_Termino2 = 0

  # Sumatoria de j a n para toda j != i
  for (j in 1:n){
    if (j != i ){

      # Operación de la sumatoria
      nbr_Termino1 = (-(mtrx_A[i,j] * vct_X[j]) / (mtrx_A[i,i]))

      # Acumulamos los valores
      nbr_Sumatoria = nbr_Sumatoria + nbr_Termino1

    }
  }

  # Terminada la sumatoria, se prepara un término extra
  nbr_Termino2 = (vct_B[i] / mtrx_A[i,i])

  # El resultado final, es lo acumulado de la sumatoria más el otro término
  nbr_Final = nbr_Sumatoria + nbr_Termino2

  # Regresamos el resultado
  nbr_Final

}

# Función que obtiene la aproximación de las iteraciones
funcObtenerVctRslt <- function(nbr_MaxIteraciones, n, mtrx_A, vct_B, vct_X0){

  # Inicializamos los vectores de control
  vct_X_Act <- vct_X0
  vct_X_Ant <- vct_X0

  # Los siguientes prints son para debuguear, más adelante se eliminarán
  print(paste0('Iteracion ', 0))
  print(vct_X_Act)

  # Máximo número de iteraciones
  for (it in 1:nbr_MaxIteraciones){

    # Iteraciones para obtener cada componente del vector de resultados
    for (i in 1:n){
      vct_X_Act[i]=funcObtenerComponente(i, n, mtrx_A, vct_X_Ant, vct_B)
    }

    # El vector resultado (k), lo usamos como vector anterior (k-1) para la sigueinte
    # iteraación
    vct_X_Ant <- vct_X_Act

    # Los siguientes prints son para debuguear, más adelante se eliminarán
    print(paste0('Iteracion ', it))
    print(vct_X_Act)

  }

  # Devolvemos el último vector calculado
  vct_X_Act

}

########################## Declaración de variables ##########################

# Construcción de los objetos matemáticos
# Matriz A
mtrx_A <- matrix(c(10,-1,2,0,
                   -1,11,-1,3,
                    2,-1,10,-1,
                     0,3,-1,8),
                 nrow=4,
                 ncol=4)

# Vector b
vct_B <- c(6, 25, -11, 15)

# N
n <- 4

# Vector de aproximación inicial
vct_X0 <- c(0,0,0,0)

# Máximo número de iteraciones permitido
nbr_MaxIteraciones <- 10

########################## Flujo principal ##########################

print('Matriz A:')
print(mtrx_A)

print('Vector b:')
print(vct_B)

# Se manda a llamar la función que obtendrá la aproximación
vct_XRslt <- funcObtenerVctRslt(nbr_MaxIteraciones, n, mtrx_A, vct_B, vct_X0)

# Se imprime el resultado
print('Resultado final: ')
print(vct_XRslt)
