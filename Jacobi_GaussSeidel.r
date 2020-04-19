########################## Información general ##########################
# Los métodos de Jacobi y Gauss_Seidel, permiten resolver
# sistemas de ecuaciones del tipo Ax=b.
# Esto se logra mediante un proceso iterativo.
# El proceso va a iterar hasta que la diferencia entre el vector de resultados
# de la iteración actual y el vector de resultados de la iteración anterior,
# sea menor a un threshold dado o que se alcance el tope de iteraciones.

# Instalación de librerías y paqueterías auxiliares
if(!require(pracma)){
    install.packages('pracma', repos = "http://cran.us.r-project.org")
}
if(!require(matrixcalc)){
    install.packages("matrixcalc")
}
library('matrixcalc')
library('pracma')
library('Matrix')

########################## Declaración de funciones ##########################

funcEsVectorValido <- function(mtrx, vct){
  # Valida si el vector (ya sea de resultados o aproximaciones) cuenta
  # con la misma cantidad de filas que la matriz a evaluar.
  #
  # Parámetros
  # ----------
  # mtrx : matriz
  #    La matriz a evaluar
  # vct : vector
  #    El vector cuyo número de filas se comparará contra las de la matriz a evaluar
  #
  # Regresa
  # -------
  # boolean
  #    Una bandera a manera de booleano indicando si la cantidad de filas del
  #    vector es igual o no a la cantidad de filas de la matriz a evaluar.
  #

  # El vector b debe tener la misma cantidad de filas que la matriz
  bool_VectorValido = FALSE
  if (nrow(mtrx) == length(vct)){
    bool_VectorValido = TRUE
  }

  bool_VectorValido

}

# Valida que sea cuadrada la matriz
funcEsMatrizCuadrada <- function(mtrx){
  # Valida si la matriz es cuadrada (nxn)
  #
  # Parámetros
  # ----------
  # mtrx : matriz
  #    La matriz a evaluar
  #
  # Regresa
  # -------
  # boolean
  #    Una bandera a manera de booleano indicando si la cantidad de filas y
  #    columnas de la matriz son iguales.
  #

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
  # Valida si se presenta algún cero en la diagonal de la matriz
  #
  # Parámetros
  # ----------
  # mtrx : matriz
  #    La matriz a evaluar
  #
  # Regresa
  # -------
  # boolean
  #    Una bandera a manera de booleano indicando si se encontró un cero en
  #    la diagonal de la matriz.
  #

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
  # Obtiene un solo componente del vector de aproximaciones
  #
  # Parámetros
  # ----------
  # i : número
  #    Indíce del componente (del vector de aproximaciones) que se quiere obtener
  # n : número
  #    Dimensión de filas o columnas de la matriz
  # mtrx_A: matriz
  #    La matriz que se está evaluando
  # vct_X : vector
  #    Vector de aproximaciones
  # vct_B : vector
  #    Vector de resultados del sistema de ecuaciones
  #
  # Regresa
  # -------
  # número
  #    El valor del componente indicado para el vector de aproximaciones
  #

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

  # El resultado final es lo acumulado de la sumatoria más el otro término
  nbr_Final = nbr_Sumatoria + nbr_Termino2

  # Regresamos el resultado
  nbr_Final

}

# Función que obtiene la aproximación de las iteraciones
funcObtenerVctRslt <- function(nbr_MaxIteraciones, n, mtrx_A, vct_B, vct_X0, nbr_Threshold, str_Metodo){
  # Mediante un proceso de iteraciones, actualiza el vector de aproximaciones
  # vct_X0 para lograr la igualdad: mtrx_A * vct_X0 = vct_B
  # El proceso de iteraciones está sujeto a que se cumpla alguna de las
  # siguientes 2 condiciones:
  #   1: Alcanzar el máximo número de iteraciones (especificado en el parámetro
  #     nbr_MaxIteraciones)
  #   2: Lograr llegar a un diferencia entre iteraciones menor al threshold
  #     que se especifica en el parámetro nbr_Threshold.
  # En cuanto se cumpla alguna de dichas condiciones, termina el proceso de
  # iteraciones.
  # Adicionalmente, hay 2 maneras de calcular el vector de resultados, mediante
  # el método Jacobi o Gauss-Seidel. La manera de especificar qué método se
  # quiere emplear es con el parámetro: str_Metodo que se emplea de la
  # siguiente manera:
  # Si el valor de str_Metodo es 'J', se emplea el método Jacobi
  # Si el valor de str_Metodo es 'GS', se emplea el método Gauss-Seidel
  #
  # Parámetros
  # ----------
  # nbr_MaxIteraciones : número
  #    Máximo número de iteraciones a alcanzar
  # n : Número
  #    Dimensión de filas o columnas de la matriz
  # mtrx_A : matriz
  #    Matriz a evaluar
  # vct_B : vector
  #    Vector de resultados del sistema de ecuaciones
  # vct_X0 : vector
  #    Vector de aproximaciones
  # nbr_Threshold : número
  #    Diferencia mínima a la que se quiere llegar entre iteraciones para
  #    considerar que el método ha convergido
  # str_Metodo : cadena
  #    Cadena mediante la cual se especifica el método que se empleará para
  #    actualizar el vector de aproximaciones.
  #
  # Regresa
  # -------
  # vector
  #    El vector de aproximaciones actualizado luego del proceso de iteraciones
  #

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

    # Si se empiezan a obtener valores NaN, significa que no hay solución
    if (is.na(nbr_Diff)){
      print('No hay convergencia')
      vct_X_Act <- rep(NA, size(vct_X0)[2])
      break
    }

    # Si se llega a una diferencia menor al threshold indicado, salimos del for
    if (nbr_Diff<nbr_Threshold){
      print('Se alcanza el threshold')
      break
    }

    # El vector resultado (k) lo usamos como vector anterior (k-1) para la sigueinte
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
  # Intercambia las filas de un vector
  #
  # Parámetros
  # ----------
  # mtrx_A : vector
  #    El vector donde se realizará el intercambio de filas
  # nbr_FilaOrigen : número
  #    Índice de la fila origen que se moverá a la fila destino
  # nbr_FilaDestino : número
  #    Índice de la fila destino que se moverá a la fila origen
  #
  # Regresa
  # -------
  # vector
  #    El vector con los valores intercambiados
  #

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
  # Intercambia las filas de una matriz
  #
  # Parámetros
  # ----------
  # mtrx : matriz
  #    La matriz donde se realizará el intercambio de filas
  # nbr_FilaOrigen : número
  #    Índice de la fila origen que se moverá a la fila destino
  # nbr_FilaDestino : número
  #    Índice de la fila destino que se moverá a la fila origen
  #
  # Regresa
  # -------
  # matriz
  #    La matriz con los valores intercambiados
  #

  # print('funcInterCambiarFilasMtrx')
  # print(nbr_FilaOrigen)
  # print(nbr_FilaDestino)

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
  # Ordena las ecuaciones del sistema buscando que no quede ningún cero
  # sobre la diagonal principal (no hay garantía de que no quede algún
  # cero sobre la diagonal). Para saber qué fila tomar, se hace una búsqueda
  # sobre cada columna preguntando por la norma infinito de cada vector-columna.
  # Puesto que el sistema de ecuaciones consta tanto de variables como de
  # resultados, es necesario tambén el re-acomodo del vector de resultados.
  #
  # Parámetros
  # ----------
  # mtrx_A : matriz
  #    La matriz que se va a ordenar
  # vct_B : vector
  #    El vector que se va a ordenar.
  #
  # Regresa
  # -------
  # list
  #    Una lista que contiene la matriz a evaluar y el vector de resultados.
  #

  # Variables que se usan dentro de la función
  nbr_Filas <- nrow(mtrx_A)
  nbr_Cols <- ncol(mtrx_A)

  bool_NrmInfPos <- TRUE

  # Se barren todas las columnas (iterador j)
  for (j in 1:(nbr_Cols-1)){
    #print('Inicia iteracion')
    #print(j)
    bool_NrmInfPos <- TRUE

    # Se saca el vector-columna que se usará en esta iteración
    vct_Col = mtrx_A[j:nbr_Cols,j]

    # Mostramos el vector-columna con el que trabajaremos
    # print(vct_Col)

    # Se obtiene la norma infinita del vector-columna
    nbr_Norm = Norm(vct_Col, p = Inf)
    # print(nbr_Norm)

    # Se pregunta si el valor es único en el vector-columna
    vct_OrdDesc <- sort(vct_Col, decreasing = TRUE)
    # print(vct_OrdDesc)

    # Preguntamos si el primer elemento del vector es diferente a la norma infinito
    # (eso significará que la norma infinito proviene de un número negativo)
    if (nbr_Norm != vct_OrdDesc[1]){
      vct_OrdDesc <- sort(vct_Col, decreasing = FALSE)
      bool_NrmInfPos <- FALSE
    }

    # Si sí es único:
    if (vct_OrdDesc[1] != vct_OrdDesc[2]){
      #print('Es unico')

      # Se obtiene el índice donde está ese valor
      nbr_Index <- match(nbr_Norm,vct_Col)
      # print(nbr_Index)

      # Si el índice es un NA
      if (is.na(nbr_Index)==TRUE){

        # Multiplicamos el valor de la norma infinito por -1
        nbr_Index <- match(nbr_Norm * -1,vct_Col)

      }

      nbr_Index <- nbr_Index + (j-1)
      # Para realizar el intercambio de filas, nuestra
      # fila origen será: nbr_Index, y la fila destino: j
      mtrx_A <- funcInterCambiarFilasMtrx(mtrx_A, nbr_Index, j, nbr_Cols)
      vct_B <- funcInterCambiarFilasVct(vct_B, nbr_Index, j)

      # print(mtrx_A)

    } else { # Si no es único:

      # Preguntamos si la norma infinito viene de un positivo o negativo
      if (bool_NrmInfPos==FALSE){
        # La multiplicamos por -1 para que haya coincidencias
        nbr_Norm = nbr_Norm * -1
      }

      #if (j==94) {
        #print('Hay empate')

        #print('mtrx_A')
        #print(mtrx_A)

        #print('nbr_Norm')
        #print(nbr_Norm)

        #print('vct_Col')
        #print(vct_Col)
      #}

      vct_Bool1 <- (vct_Col==nbr_Norm)

      # Si se trata de la última columna, es un caso especial
      if ((j+1)==nbr_Cols){
        #print('Caso especial')
        mtrx_Tmp1 <- mtrx_A[j:nbr_Cols, (j+1):nbr_Cols]
        mtrx_Tmp1 <- rbind(mtrx_Tmp1)
        mtrx_Tmp1 <- t(mtrx_Tmp1)

        mtrx_Tmp2 <- mtrx_Tmp1[vct_Bool1, 1:(nbr_Cols-j)]
        mtrx_Tmp2 <- rbind(mtrx_Tmp2)
        mtrx_Tmp2 <- t(mtrx_Tmp2)

      }else{ # Este es el caso normal
        mtrx_Tmp1 <- mtrx_A[j:nbr_Cols, (j+1):nbr_Cols]
        mtrx_Tmp2 <- mtrx_Tmp1[vct_Bool1, 1:(nbr_Cols-j)]
      }

      #if (j==94) {
        #print('mtrx_Tmp1 pt1')
        #print(mtrx_Tmp1)
        #print('dim')
        #print(dim(mtrx_Tmp1))
        #print('vct_Bool1')
        #print(vct_Bool1)
        #print('mtrx_Tmp2 pt2')
        #print(mtrx_Tmp2)
        #print('dim')
        #print(dim(mtrx_Tmp2))
      #}


      # Barremos el resto de las columnas para el desempate
      for (jj in 1:(nbr_Cols-j)){
        #print('j')
        #print(j)
        #print('jj')
        #print(jj)
        #print('nrow(mtrx_Tmp2)')
        #print(nrow(mtrx_Tmp2))
        vct_ColDesempate <- mtrx_Tmp2[jj:nrow(mtrx_Tmp2),jj]

        # Buscaremos el valor mínimo de la siguiente columna
        vct_OrdAsc <- sort(vct_ColDesempate)
        #print(vct_OrdAsc)

        # Si el  mínimo es único:
        if (vct_OrdAsc[1] != vct_OrdAsc[2]){

          #print('es unico')
          # Generamos el vector que nos ayudará a obtener las filas a desempatar
          vct_Bool2 <- (vct_ColDesempate == vct_OrdAsc[1])
          #print(vct_Bool2)

          # Ya que tenemos el vector que dice cuál es la fila que debemos tomar
          # vamos yendo hacia atrás para averiguar el índice relativo de
          # en columna con la que estamos trabajando
          #print('Comienza regreso')
          vct_Tmp1 <- mtrx_Tmp2[vct_Bool2]
          #print('vct_Tmp1')
          #print(vct_Tmp1)

          vct_Bool3 <- (vct_Tmp1==mtrx_Tmp1)
          #print('vct_Bool3')
          #print(vct_Bool3)

          vct_bool4 <- apply(mtrx_Tmp1, 1, function(x) all(x==vct_Tmp1))
          #print('vct_bool4')
          #print(vct_bool4)

          nbr_Index <- match(TRUE, vct_bool4)
          # nbr_Index <- match(TRUE,vct_Bool2)
          #print('nbr_Index')
          #print(nbr_Index)
          nbr_Index <- nbr_Index + (j-1)

          # Para realizar el intercambio de filas, nuestra
          # fila origen será: nbr_Index y la fila destino: j
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
  # Resuelve un sistema de ecuaciones lineales mediante el método de Jacobi
  # o de Gauss-Seidel. El sistema de ecuaciones sólo se procesará si pasa
  # todas las validaciones requeridas.
  #
  # Parámetros
  # -------
  # mtrx_A : matriz
  #    Matriz a evaluar
  # vct_B : vector
  #    Vector de resultados del sistema de ecuaciones
  # vct_X0 : vector
  #    Vector de aproximaciones
  # nbr_MaxIteraciones : número
  #    Máximo número de iteraciones a alcanzar
  # nbr_Threshold : número
  #    Diferencia mínima a la que se quiere llegar entre iteraciones para
  #    considerar que el método ha convergido
  # str_Metodo : cadena
  #    Cadena mediante la cual se especifica el método que se empleará para
  #    actualizar el vector de aproximaciones.
  #
  # Regresa
  # -------
  # vector
  #    El vector de aproximaciones luego del proceso de iteraciones

  # Se agrega condicón para validar que la matriz no sea singular

  # Se inicializa el vector de resultados
  vct_XRslt <- rep(NA, size(vct_X0)[2])

  if(!is.singular.matrix(mtrx_A)){

    if (str_Metodo == 'J' || str_Metodo == 'GS'){

      if (str_Metodo=='J'){
        print('Solución mediante metodo de Jacobi')
      }
      if (str_Metodo=='GS'){
        print('Solución mediante metodo de Gauss-Seidel')
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
              print('La matriz tiene algún cero en la diagonal, comienza ordenamiento')

              # Respaldo del orden original
              mtrx_Original <- mtrx_A
              vct_Original <- vct_B

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
                print('Resultado final (nuevo orden): ')
                print(vct_XRslt)

                # Se reordena el vector de resultados para entregarlo en el mismo
                # orden en que fue solicitado
                nbr_Filas <-size(vct_XRslt)[2]
                for (i in 1:nbr_Filas){
                  for (j in 1:nbr_Filas){
                    if (sum(mtrx_A[i,]==mtrx_Original[j,]) == nbr_Filas){
                      mtrx_A <- funcInterCambiarFilasMtrx(mtrx_A, i, j, nbr_Filas)
                      vct_XRslt <- funcInterCambiarFilasVct(vct_XRslt, i, j)
                    }
                  }
                }
                print('Resultado final (orden original): ')
                print(vct_XRslt)

              } else {
                print('Pese al reordenamiento, aún hay ceros en la diagonal')
              }

            }

          } else {
              print('La matriz no cumple con ser de dimensiones nxn')
          }
        } else {
          print('El vector de aproximaciones no es de las dimensiones esperadas')
        }
      } else {
        print('El vector de resultados no es de las dimensiones esperadas')
      }

    }else{
      print("El método especificado no es válido, se espera 'GS' para Gauss-Seidel o 'J' para Jacobi")
    }
  }else{
    print('La matriz es singular, por lo tanto no hay solución al sistema y el método se detiene')
  }

  # Se devuelve el vector de resultados
  vct_XRslt

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
# vct_B <- c(6, 25, -11, 15)

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
# del método de Jacobi o Gauss-Seide, según se indique
vct_Solucion <- funcResolverSE(mtrx_A, vct_B, vct_X0, nbr_MaxIteraciones, nbr_Threshold, str_Metodo)
