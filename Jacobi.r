########################## Información general ##########################
# Los métodos de Jacobi y Gauss_Seidel, permiten resolver
# sistemas de ecuaciones del tipo Ax=b.
# Esto se logra mediante un proceso iterativo.
# El proceso va a iterar hasta que la diferencia entre el vector de resultados
# de la iteración actual y el vector de resultados de la iteración anterior,
# sea menor a un threshold dado o que se alcance el tope de iteraciones.

# Instalación de librerías y paqueterías auxiliares
install.packages('pracma', repos = "http://cran.us.r-project.org")
library('pracma')

########################## Declaración de funciones ##########################

funcEsVectorValido <- function(mtrx, vct){

  # El vector b, debe tener la misma cantidad de filas que la matriz
  bool_VectorValido = FALSE
  if (nrow(mtrx) == length(vct)){
    bool_VectorValido = TRUE
  }

  bool_VectorValido

}

# Valida que sea cuadrada la matriz
funcEsMatrizCuadrada <- function(mtrx){

  # Si el número de filas es igual al número de columnas, es una matriz cuadrada
  if (nrow(mtrx)==ncol(mtrx)){
    bool_Valida <- TRUE
  }
  else{
    bool_Valida <- FALSE
  }

  bool_Valida

}

 # Busca algún cero en la diagonal principal
funcHayCeroEnDiagonal <- function(mtrx){

  bool_HayCero <- FALSE
  nbr_Filas <- nrow(mtrx)
  nbr_Cols <- ncol(mtrx)

  for (i in 1:nbr_Filas){
    for (j in 1:nbr_Cols){
      if ((i==j) && (mtrx[i,j]==0)){
        bool_HayCero <- TRUE
      }
    }
  }

  bool_HayCero

}

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
funcObtenerVctRslt <- function(nbr_MaxIteraciones, n, mtrx_A, vct_B, vct_X0, nbr_Threshold, str_Metodo){

  # Inicializamos los vectores de control
  vct_X_Act <- vct_X0
  vct_X_Ant <- vct_X0

  # Los siguientes prints son para debuguear, más adelante se eliminarán
  print(paste0('Iteracion ', 0))
  print(vct_X_Act)

  # Máximo número de iteraciones
  for (it in 1:nbr_MaxIteraciones){

    print(paste0('Iteracion ', it))

    # Iteraciones para obtener cada componente del vector de resultados
    for (i in 1:n){

      # Si se pidió usar el método Jacobi
      if (str_Metodo=='J'){
        vct_X_Act[i]=funcObtenerComponente(i, n, mtrx_A, vct_X_Ant, vct_B)
      }

      # Si se pidió usar el método Gauss-Seidel
      if (str_Metodo=='GS'){
        vct_X_Act[i]=funcObtenerComponente(i, n, mtrx_A, vct_X_Act, vct_B)
      }

    }

    print(vct_X_Act)

    nbr_Numerador <- Norm(vct_X_Act - vct_X_Ant, p = Inf)
    nbr_Denominador <- Norm(vct_X_Act, p = Inf)

    print(paste0('nbr_Numerador: ', nbr_Numerador))
    print(paste0('nbr_Denominador: ', nbr_Denominador))

    nbr_Diff <-  nbr_Numerador / nbr_Denominador
    print(paste0('nbr_Diff: ',nbr_Diff))

    # Si se llega a una diferencia menor al threshold indicado, salimos del for
    if (nbr_Diff<nbr_Threshold){
      print('Se alcanza el threshold')
      break
    }

    # El vector resultado (k), lo usamos como vector anterior (k-1) para la sigueinte
    # iteraación
    vct_X_Ant <- vct_X_Act

  }

  if (it==nbr_MaxIteraciones){
    print('Se llega al tope de iteraciones')
  }

  # Devolvemos el último vector calculado
  vct_X_Act

}

funcInterCambiarFilasVct <- function(vctOrigen, nbr_FilaOrigen, nbr_FilaDestino){

  # Se guarda el valor destino
  nbr_ValorTmp <- vctOrigen[nbr_FilaDestino]

  # Se pone el valor origen hacia el valor destino
  vctOrigen[nbr_FilaDestino] <- vctOrigen[nbr_FilaOrigen]

  # Se recupera el valor destino original, y se pone en valor origen
  vctOrigen[nbr_FilaOrigen] <- nbr_ValorTmp

  # Se regresa el valor
  vctOrigen

}

funcInterCambiarFilasMtrx <- function(mtrx, nbr_FilaOrigen, nbr_FilaDestino, nbr_Cols){

  # Se guarda el valor destino
  vct_FilaTmp <- mtrx[nbr_FilaDestino,1:nbr_Cols]

  # Se pone el valor origen hacia el valor destino
  mtrx[nbr_FilaDestino,1:nbr_Cols] <- mtrx[nbr_FilaOrigen,1:nbr_Cols]

  # Se recupera el valor destino original, y se pone en valor origen
  mtrx[nbr_FilaOrigen,1:nbr_Cols] <- vct_FilaTmp

  # Se regresa el valor
  mtrx

}

funcOrdenarEcuaciones <- function(mtrx_A, vct_B){

  # Variables que se usan dentro de la funcion
  nbr_Filas <- nrow(mtrx_A)
  nbr_Cols <- ncol(mtrx_A)

  # Se barren todas las columnas (iterador j)
  for (j in 1:nbr_Cols){
    # print('Inicia iteracion')

    # Se saca el vector-columna que se usará esta iteración
    vct_Col = mtrx_A[1:nbr_Cols,j]

    # Mostramos el vector-columna con el que trabajaremos
    # print(vct_Col)

    # Se obtiene la norma infinita del vector-columna
    nbr_Norm = Norm(vct_Col, p = Inf)
    # print(nbr_Norm)

    # Se pregunta si el valor es único en el vector-columna
    vct_OrdDesc <- sort(vct_Col, decreasing = TRUE)
    # print(vct_OrdDesc)

    # Si sí es único:
    if (vct_OrdDesc[1] != vct_OrdDesc[2]){
      #print('Es unico')
      # Se obtiene el índice donde está ese valor
      nbr_Index <- match(nbr_Norm,vct_Col)
      # print(nbr_Index)

      # Para realizar el intercambio de filas, nuestra
      # fila origen será: nbr_Index, y la fila destino: j
      mtrx_A <- funcInterCambiarFilasMtrx(mtrx_A, nbr_Index, j, nbr_Cols)
      vct_B <- funcInterCambiarFilasVct(vct_B, nbr_Index, j)

      # print(mtrx_A)

    } else { # Si no es único:
      #print('Hay empate')

      vct_Bool1 <- (vct_Col==nbr_Norm)

      mtrx_Tmp <- mtrx_A[vct_Bool1,1:nbr_Cols]
      #print(mtrx_Tmp)

      # Barremos el resto de las columans para el desempate
      for (jj in (j+1):nbr_Cols){

        vct_ColDesempate <- mtrx_Tmp[,jj]
        # print(paste0('vct_ColDesempate: ', vct_ColDesempate))

        # Buscaremos el valor mínimo de la siguiente columna
        vct_OrdAsc <- sort(vct_ColDesempate)
        # print(vct_OrdAsc)

        # Si el  mínimo es único:
        if (vct_OrdAsc[1] != vct_OrdAsc[2]){

          # print('es unico')
          # Generamos el vector que nos ayudará a obtener las filas a desempatar
          vct_Bool2 <- (vct_ColDesempate == vct_OrdAsc[1])
          # print(vct_Bool2)

          vct_FilaDesempate <- mtrx_Tmp[vct_Bool2]

          # Ya que se tiene la fila mínima del empate, se obtiene su indice
          # en la matriz_A
          vct_Bool3 <- rowSums(mtrx_A == vct_FilaDesempate[col(mtrx_A)]) == ncol(mtrx_A)
          nbr_Index <- match(TRUE,vct_Bool3)

          # Para realizar el intercambio de filas, nuestra
          # fila origen será: nbr_Index, y la fila destino: j
          mtrx_A <- funcInterCambiarFilasMtrx(mtrx_A, nbr_Index, j, nbr_Cols)
          vct_B <- funcInterCambiarFilasVct(vct_B, nbr_Index, j)

          #print(mtrx_A)
          break

        }

      }

    }

  }

  # Regresamos en una lista, la matriz y vector ordenados
  list(matriz=mtrx_A,vector=vct_B)

}

funcResolverSE <- function(mtrx_A, vct_B, vct_X0, nbr_MaxIteraciones, nbr_Threshold, str_Metodo){

  if (str_Metodo == 'J' || str_Metodo == 'GS'){

    if (str_Metodo=='J'){
      print('Solucion mediante metodo de Jacobi')
    }
    if (str_Metodo=='GS'){
      print('Solucion mediante metodo de Gauss-Sidel')
    }

    print('Matriz A:')
    print(mtrx_A)

    print('Vector b:')
    print(vct_B)

    # Se aplican las validaciones de manera anidada
    if (funcEsVectorValido(mtrx_A, vct_B)){

      if (funcEsVectorValido(mtrx_A, vct_X0)){

        if (funcEsMatrizCuadrada(mtrx_A) == TRUE){

          # Obtenemos la n de la matriz
          n <- nrow(mtrx_A)

          if (funcHayCeroEnDiagonal(mtrx_A) == FALSE){

            # Se manda a llamar la función que obtendrá la aproximación
            vct_XRslt <- funcObtenerVctRslt(nbr_MaxIteraciones, n, mtrx_A, vct_B, vct_X0, nbr_Threshold, str_Metodo)

            # Se imprime el resultado
            print('Resultado final: ')
            print(vct_XRslt)

          } else {
            print('La matriz tiene algun cero en la diagonal, comienza ordenamiento')

            lt_Obj <- funcOrdenarEcuaciones(mtrx_A, vct_B)
            mtrx_A <- lt_Obj$matriz
            vct_B <- lt_Obj$vector

            print('Matriz ordenada:')
            print(mtrx_A)

            print('Vector ordenado:')
            print(vct_B)

            if (funcHayCeroEnDiagonal(mtrx_A) == FALSE){

              # Se manda a llamar la función que obtendrá la aproximación
              vct_XRslt <- funcObtenerVctRslt(nbr_MaxIteraciones, n, mtrx_A, vct_B, vct_X0, nbr_Threshold, str_Metodo)

              # Se imprime el resultado
              print('Resultado final: ')
              print(vct_XRslt)

            } else {
              print('Pese al reordenamiento, aun hay ceros en la diagonal')
            }

          }

        } else {
            print('La matriz no cumple con ser de dimensiones nxn')
        }
      } else {
        print('El vector de aproximaciones, debe tener la misma cantidad de filas que la matriz a evaluar')
      }
    } else {
      print('El vector de resultados, debe tener la misma cantidad de filas que la matriz a evaluar')
    }

  }else{
    print('El metodo especificado no es valido, favor de verificar')
  }
}


########################## Declaración de variables ##########################

# Construcción de los objetos matemáticos

# Ejemplo de Erick (en desorden)
# Matriz A
mtrx_A <- matrix(c( 0,3,-1,8,
                   -1,11,-1,3,
                    2,-1,10,-1,
                    10,-1,2,0),
                 byrow = TRUE,
                 nrow=4,
                 ncol=4)

vct_B <- c(15, 25, -11, 6)
#vct_B <- c(6, 25, -11, 15)

# Ejemplo que puso Ana en el chat
# Matriz A
# mtrx_A <- matrix(c(1,0,0,0,
#                   0,1,2,1,
#                   1,2,0,0,
#                   1,1,1,0),
#                   nrow=4,
#                   ncol=4,
#                  byrow = TRUE)

# Vector b
# vct_B <- c(10, 24, 31, 45)

# Vector de aproximación inicial
vct_X0 <- c(0,0,0,0)

# Máximo número de iteraciones permitido
nbr_MaxIteraciones <- 100

# Threshold que se busca alcanzar
nbr_Threshold <- 10**(-3)

# Método a utilizar
#str_Metodo <- 'J'
str_Metodo <- 'GS'

########################## Flujo principal ##########################

# Se manda a llamar la función que contiene las validaciones y ejecución
# del método de Jacobi
funcResolverSE(mtrx_A, vct_B, vct_X0, nbr_MaxIteraciones, nbr_Threshold, str_Metodo)
