
# Examen de cómputo matricial equipo Jacobi, Gauss-Seidel


## Equipo de trabajo

El equipo está integrado por seis personas: un project manager, dos programadores y tres revisores. 

| Alumno | Equipo |
|--------|--------|
| Marco  | Programador |
| Mario  | Programador |
| Elizabeth | Revisor  |
| Itzel | Revisor |
| Oscar | Revisor  |
| Ana   | Project Manager |

## Metodología de trabajo

El proyecto se desarrolló bajo una metodología ágil siguiendo el marco de trabajo de [scrum](https://github.com/mno-2020-gh-classroom/ex-modulo-3-comp-matricial-jacobi-gauss-seidel-anabco/blob/master/Intro-scrum.md).

Se dividió en cuatro milestones: 

+ Definición de equipo y forma de trabajo. El objetivo es definir los roles y responsabilidades. El project manager debe explicar cómo funciona scrum así como el project board. 
+ Milestone1: método de Jacobi: entendimiento, programación y pruebas del método. 
+ Milestone2: método de Gauss Seidel: entendimiento, programación y pruebas del método. 
+ Milestone3: método de eliminación por bloques: entendimiento, programación y pruebas del método. 
+ Milestone4: Documentación: documentación de los 3 métodos y revisiones finales de los archivos. 

Las tareas e issues detectados durante el desarrollo del proyecto están organizadas en un [project board](https://github.com/mno-2020-gh-classroom/ex-modulo-3-comp-matricial-jacobi-gauss-seidel-anabco/projects/1)

## Estructura del repositorio

**Archivos**
+ Parámetros_modelos: es unan nota ráoida que indica los parámetros que los programadores deberían considerar en el desarrollo de los programas. 


**Carpetas**
+ Revisiones: se realizó un documento por cada método. En la carpeta están disponibles los archivos en formato .md y .htlm. Sin embargo, para la facilitar la lectura de los documentos desde GitHub, se pueden consultar los siguientes links: 

+ [Jacobi](https://mno-2020-gh-classroom.github.io/ex-modulo-3-comp-matricial-jacobi-gauss-seidel-anabco/Revisiones/pruebas_jacobi.html)
+ [Gauss Seidel](https://mno-2020-gh-classroom.github.io/ex-modulo-3-comp-matricial-jacobi-gauss-seidel-anabco/Revisiones/pruebas_gauss_seidel.html)
+ [Eliminación por bloques]

Se crea una liga de Binder para realizas las pruebas de los 3 métodos. 

**Liga de Binder para las pruebas del programa de Jacobi, Gauss-Seidel y Eliminación por Bloques:** [![Binder](https://mybinder.org/badge_logo.svg)](https://mybinder.org/v2/gh/shimanteko/for_binders/master?urlpath=lab/tree/home/jovyan/Jacobi_GaussSeidel_Bloques.ipynb)

## Objetivos:

* Programar el método de eliminación por bloques que se encuentra en la sección **Métodos o algoritmos numéricos por bloques para SEL** en la nota [3.3.Solucion_de_SEL_y_FM](https://github.com/ITAM-DS/analisis-numerico-computo-cientifico/blob/master/temas/III.computo_matricial/3.3.Solucion_de_SEL_y_FM.ipynb). Para resolver este método, el equipo también programará [3.3.e.Jacobi_Gauss-Seidel](https://github.com/ITAM-DS/analisis-numerico-computo-cientifico/blob/master/temas/III.computo_matricial/3.3.e.Jacobi_Gauss-Seidel.ipynb) para resolver los sistemas de ecuaciones que surjan en el método de eliminación por bloques. Utilizar matrices pseudoaleatorias de tamaño mediano: aprox de dimensiones de $10^4 \times 10^4$.

* Aprendizaje sobre el uso de github como herramienta colaborativa en la creación y desarrollo de proyectos.

* Aprendizaje en la organización de trabajo en equipo para adopción de frameworks como [scrum](https://www.youtube.com/watch?v=b02ZkndLk1Y&feature=emb_logo) para el desarrollo de proyectos. 

**Nota para los equipos que programan una factorización matricial distinta a la SVD:** una vez que programan su factorización matricial tienen que utilizar métodos de sustitución hacia delante y hacia atrás para resolver el SEL asociado. Utilicen [solve triangular](https://docs.scipy.org/doc/scipy/reference/generated/scipy.linalg.solve_triangular.html) de `scipy` o [backsolve/forwardsolve](https://stat.ethz.ch/R-manual/R-devel/library/base/html/backsolve.html) de `R` para esto.

## Fecha de entrega

* 19 de abril 11:59 pm

## Lenguaje a utilizar: R




