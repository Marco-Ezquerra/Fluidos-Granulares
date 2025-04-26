Simulación fluidos granulares

Este proyecto implementa una simulación de un gas granular en régimen diluido, donde las partículas sufren colisiones inelásticas modeladas mediante una versión simplificada de la ecuación de Boltzmann. Se analiza la evolución de las velocidades y los momentos cinéticos de las partículas bajo un proceso de enfriamiento homogéneo.

📋 Descripción del modelo

- Sistema de N = 131072 partículas (2¹⁷).
- Evolución en el espacio 3D (componentes de velocidad: vx, vy, vz).
- Colisiones binarias con coeficiente de restitución a = 0.8.
- Condiciones iniciales: distribución direccional aleatoria.
- Ajuste periódico de momento total (con reajuste) para evitar deriva numérica.
- Cálculo de momentos de velocidad de segundo y cuarto orden (μ2, μ4).
- Histograma de velocidades para análisis posterior.
- Número total de colisiones: 100,000, del orden del numero de particulas

 Estructura del código

fluidos_granulares.f95     # Código principal en Fortran 95
velocidades.txt            # Archivo de salida con componentes de velocidad

 Requisitos

- Compilador Fortran 95 (recomendado: gfortran).
- Gnuplot (opcional) para graficar histogramas.

 Compilación y ejecución

gfortran -o granular fluidos_granulares.f95
./granular

El código genera el archivo velocidades.txt, que contiene las velocidades tras las colisiones, listo para ser analizado estadísticamente o visualizado.

 Observables calculados

- Momentos cinéticos:
  μ2 = <v^2>,  μ4 = <v^4>
- Histogramas de distribución de velocidades en cada dirección (hx, hy, hz).
- Velocidad típica: v_típica ~ sqrt(μ2)
- Coeficientes de cumulantes para analizar desviaciones de la distribución de Maxwell-Boltzmann.

Fundamento teórico

Este modelo sigue el marco cinético para gases granulares tal como se describe en:

- Haff, P. K. (1983). "Grain flow as a fluid-mechanical phenomenon".
- Brilliantov & Pöschel (2004). "Kinetic Theory of Granular Gases".
- Goldhirsch (2003). "Rapid granular flows", Annu. Rev. Fluid Mech.
