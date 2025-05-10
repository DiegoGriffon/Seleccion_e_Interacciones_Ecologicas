# Script en R para simular el modelo ampliado de Leon (1974)
# Selección en contextos de interacciones entre especies (competencia, mutualismo, etc.)


# Definir parámetros base para Species A (A1A1, A1A2, A2A2)
RA_params <- c(0.1, 0.1, 0.1)
KA_params <- c(100, 150, 120)
alpha_base_params <- c(1.0, 0.5, 0.8) # Valores base de alpha (siempre no negativos)

# Definir parámetros base para Species B (B1B1, B1B2, B2B2)
RB_params <- c(0.15, 0.15, 0.15)
KB_params <- c(110, 160, 130)
beta_base_params <- c(0.9, 0.4, 0.7) # Valores base de beta (siempre no negativos)

# Condiciones iniciales
initial_NA <- 50
initial_NB <- 60
initial_pA1 <- 0.5
initial_pB1 <- 0.5

# Número de generaciones a simular
num_generations <- 400

# El tipo de interacción entre las especies se asigna de forma aleatoria y se indica en las salidas graficas de la simulación.

# Define la función de simulación
# Añade alpha_multiplier y beta_multiplier como argumentos
simulate_leon_model_flexible <- function(NA_init, NB_init, pA1_init, pB1_init,
                                         RA, KA, alpha, RB, KB, beta,
                                         alpha_multiplier, beta_multiplier, # Nuevos argumentos para el tipo de interacción
                                         ngen, interaction_label = "Interacción desconocida") { # Para usar en los gráficos
  
  # --- Validar entradas básicas ---
  if (length(RA) != 3 || length(KA) != 3 || length(alpha) != 3) {
    stop("RA, KA y alpha deben ser vectores de longitud 3 (para genotipos 11, 12, 22).")
  }
  if (length(RB) != 3 || length(KB) != 3 || length(beta) != 3) {
    stop("RB, KB y beta deben ser vectores de longitud 3 (para genotipos 11, 12, 22).")
  }
  if (any(RA <= 0) || any(RA >= 1) || any(RB <= 0) || any(RB >= 1)) {
    warning("Se recomienda que los valores de R estén entre 0 y 1 para evitar oscilaciones en el modelo logístico discreto base.")
  }
  if (any(KA <= 0) || any(KB <= 0)) {
    stop("Los valores de K deben ser positivos.")
  }
  # Alpha y Beta originales pueden ser 0, pero los multiplicadores deben ser 1, 0 o -1
  if (!all(alpha >= 0) || !all(beta >= 0)) {
    warning("Se espera que los valores base de alpha y beta sean no negativos, los multiplicadores manejan los signos.")
  }
  if (!alpha_multiplier %in% c(-1, 0, 1) || !beta_multiplier %in% c(-1, 0, 1)) {
    stop("Los multiplicadores de interacción deben ser -1, 0 o 1.")
  }
  if (NA_init < 0 || NB_init < 0 || pA1_init < 0 || pA1_init > 1 || pB1_init < 0 || pB1_init > 1) {
    stop("Valores iniciales inválidos. Las poblaciones deben ser no negativas y las frecuencias génicas entre 0 y 1.")
  }
  if (NA_init == 0 && NB_init == 0) {
    stop("Las poblaciones iniciales no pueden ser ambas cero.")
  }
  
  
  # --- Inicializar almacenamiento ---
  results <- data.frame(
    generation = 0:ngen,
    NA_pop = numeric(ngen + 1),
    NB_pop = numeric(ngen + 1),
    pA1 = numeric(ngen + 1),
    pB1 = numeric(ngen + 1)
  )
  
  # Almacenar condiciones iniciales
  results$NA_pop[1] <- NA_init
  results$NB_pop[1] <- NB_init
  results$pA1[1] <- pA1_init
  results$pB1[1] <- pB1_init
  
  # --- Bucle de simulación ---
  for (t in 1:ngen) {
    # Obtener variables de estado actuales
    NA_t <- results$NA_pop[t]
    NB_t <- results$NB_pop[t]
    pA1_t <- results$pA1[t]
    pB1_t <- results$pB1[t]
    
    # Si ambas poblaciones están extintas, detener la simulación
    if (NA_t < 1e-6 && NB_t < 1e-6) {
      warning(paste("Ambas poblaciones se han extinguido en la generación", t-1, ". Deteniendo simulación."))
      results[(t+1):(ngen+1), ] <- NA # Rellenar filas restantes con NA
      break # Salir del bucle
    }
    
    # Calcular frecuencias genotípicas (asumiendo Equilibrio Hardy-Weinberg)
    # Especie A: f(A1A1) = pA1^2, f(A1A2) = 2*pA1*(1-pA1), f(A2A2) = (1-pA1)^2
    pAA11 <- pA1_t^2
    pAA12 <- 2 * pA1_t * (1 - pA1_t)
    pAA22 <- (1 - pA1_t)^2
    
    # Especie B: f(B1B1) = pB1^2, f(B1B2) = 2*pB1*(1-pB1), f(B2B2) = (1-pB1)^2
    pBB11 <- pB1_t^2
    pBB12 <- 2 * pB1_t * (1 - pB1_t)
    pBB22 <- (1 - pB1_t)^2
    
    
    # Calcular aptitudes absolutas para la generación actual (Modificación de Eq 1 & 2)
    # Se aplica el multiplicador (+1, -1 o 0) a los términos de interacción interespecífica
    
    # Especie A genotipos (A1A1, A1A2, A2A2)
    # La interacción con NB ahora puede ser sumada, restada o cero
    W_AA11 <- 1 + (RA[1] / KA[1]) * (KA[1] - NA_t + alpha_multiplier * alpha[1] * NB_t)
    W_AA12 <- 1 + (RA[2] / KA[2]) * (KA[2] - NA_t + alpha_multiplier * alpha[2] * NB_t)
    W_AA22 <- 1 + (RA[3] / KA[3]) * (KA[3] - NA_t + alpha_multiplier * alpha[3] * NB_t)
    
    # Especie B genotipos (B1B1, B1B2, B2B2)
    # La interacción con NA ahora puede ser sumada, restada o cero
    V_BB11 <- 1 + (RB[1] / KB[1]) * (KB[1] - NB_t + beta_multiplier * beta[1] * NA_t)
    V_BB12 <- 1 + (RB[2] / KB[2]) * (KB[2] - NB_t + beta_multiplier * beta[2] * NA_t)
    V_BB22 <- 1 + (RB[3] / KB[3]) * (KB[3] - NB_t + beta_multiplier * beta[3] * NA_t)
    
    
    # Asegurar que las aptitudes sean no negativas
    W_AA11 <- max(0, W_AA11)
    W_AA12 <- max(0, W_AA12)
    W_AA22 <- max(0, W_AA22)
    V_BB11 <- max(0, V_BB11)
    V_BB12 <- max(0, V_BB12)
    V_BB22 <- max(0, V_BB22)
    
    # Calcular aptitudes promedio (Mean Fitnesses) (Eq 6a & 6b)
    W_bar <- pAA11 * W_AA11 + pAA12 * W_AA12 + pAA22 * W_AA22
    V_bar <- pBB11 * V_BB11 + pBB12 * V_BB12 + pBB22 * V_BB22
    
    # Si la aptitud promedio es cero, la población no puede crecer
    if (W_bar <= 0) W_bar = 0
    if (V_bar <= 0) V_bar = 0
    
    # Calcular aptitudes marginales (Marginal Fitnesses) (Eq 5a & 5b)
    W_marg_A1 <- pA1_t * W_AA11 + (1 - pA1_t) * W_AA12
    V_marg_B1 <- pB1_t * V_BB11 + (1 - pB1_t) * V_BB12
    
    # Actualizar frecuencias génicas para la siguiente generación (Eq 3a & 3b)
    pA1_next <- pA1_t # Mantener frecuencia si W_bar es cero o población 0
    if (W_bar > 1e-9 && NA_t > 1e-9) {
      pA1_next <- pA1_t * (W_marg_A1 / W_bar)
    }
    
    pB1_next <- pB1_t # Mantener frecuencia si V_bar es cero o población 0
    if (V_bar > 1e-9 && NB_t > 1e-9) {
      pB1_next <- pB1_t * (V_marg_B1 / V_bar)
    }
    
    
    # Asegurar que las frecuencias génicas se mantengan en [0, 1]
    pA1_next <- max(0, min(1, pA1_next))
    pB1_next <- max(0, min(1, pB1_next))
    
    # Actualizar tamaños poblacionales para la siguiente generación (Eq 4a & 4b)
    NA_next <- NA_t * W_bar
    NB_next <- NB_t * V_bar
    
    # Asegurar que los tamaños poblacionales sean no negativos
    NA_next <- max(0, NA_next)
    NB_next <- max(0, NB_next)
    
    # Almacenar resultados para la siguiente generación
    results$NA_pop[t + 1] <- NA_next
    results$NB_pop[t + 1] <- NB_next
    results$pA1[t + 1] <- pA1_next
    results$pB1[t + 1] <- pB1_next
    
    # Detener si las poblaciones se vuelven efectivamente cero
    if (NA_next < 1e-3 && NB_next < 1e-3) { # Usar un umbral pequeño para extinción
      warning(paste("Ambas poblaciones efectivamente extintas en la generación", t, ". Deteniendo simulación."))
      results[(t+2):(ngen+1), ] <- NA # Rellenar filas restantes con NA
      break # Salir del bucle
    }
  }
  
  return(results)
}

# Definir los posibles tipos de interacción y sus multiplicadores correspondientes
# Formato: Lista(c(multiplicador_alpha, multiplicador_beta), label)
interaction_types <- list(
  list(c(-1, -1), "Competencia (-/-)"),
  list(c(1, 1), "Mutualismo (+/+)"),
  list(c(1, -1), "Depredación (+/-: A gana/B pierde)"), # A es el que tiene alpha
  list(c(-1, 1), "Depredación (-/+: B gana/A pierde)"), # B es el que tiene beta
  list(c(1, 0), "Comensalismo (+/0: A gana/B neutro)"),
  list(c(0, 1), "Comensalismo (0/+: A neutro/B gana)"),
  list(c(-1, 0), "Amensalismo (-/0: A pierde/B neutro)"),
  list(c(0, -1), "Amensalismo (0/-: A neutro/B pierde)")
)

# --- Asignación aleatoria del tipo de interacción para esta ejecución ---
chosen_interaction <- sample(interaction_types, 1)[[1]] # Selecciona una lista aleatoria
alpha_mult <- chosen_interaction[[1]][1]
beta_mult <- chosen_interaction[[1]][2]
interaction_label <- chosen_interaction[[2]]

cat("Simulando interacción tipo:", interaction_label, "\n") # Imprimir el tipo de interacción

# Ejecutar la simulación con el tipo de interacción elegido
sim_results <- simulate_leon_model_flexible(
  NA_init = initial_NA,
  NB_init = initial_NB,
  pA1_init = initial_pA1,
  pB1_init = initial_pB1,
  RA = RA_params,
  KA = KA_params,
  alpha = alpha_base_params, # Pasamos los valores base positivos
  RB = RB_params,
  KB = KB_params,
  beta = beta_base_params,   # Pasamos los valores base positivos
  alpha_multiplier = alpha_mult, # Pasamos el multiplicador elegido
  beta_multiplier = beta_mult,   # Pasamos el multiplicador elegido
  ngen = num_generations,
  interaction_label = interaction_label # Pasamos la etiqueta para los gráficos
)

# --- Visualizar Resultados ---

# Librería necesaria para graficar
library(ggplot2)

# Gráfico de tamaños poblacionales
pop_plot <- ggplot(sim_results, aes(x = generation)) +
  geom_line(aes(y = NA_pop, color = "Especie A")) +
  geom_line(aes(y = NB_pop, color = "Especie B")) +
  labs(title = paste("Dinámica Poblacional:", interaction_label),
       y = "Tamaño Poblacional",
       x = "Generación",
       color = "Especie") +
  theme_minimal()

print(pop_plot)

# Gráfico de frecuencias génicas
gene_freq_plot <- ggplot(sim_results, aes(x = generation)) +
  geom_line(aes(y = pA1, color = "Especie A (pA1)")) +
  geom_line(aes(y = pB1, color = "Especie B (pB1)")) +
  labs(title = paste("Dinámica de Frecuencias Génicas:", interaction_label),
       y = "Frecuencia del Alelo 1",
       x = "Generación",
       color = "Especie/Alelo") +
  theme_minimal() +
  ylim(0, 1) # Asegurar que el eje y esté entre 0 y 1

print(gene_freq_plot)

# Gráfico de fase (Tamaño poblacional de A vs Tamaño poblacional de B)
# Solo si la simulación no terminó inmediatamente por extinción
if (sum(!is.na(sim_results$NA_pop)) > 1) {
  phase_plot <- ggplot(sim_results, aes(x = NA_pop, y = NB_pop)) +
    geom_path(color = "blue", arrow = arrow(type = "open", length = unit(0.1, "inches"))) +
    geom_point(data = sim_results[1,], color = "red", size = 2) + # Punto de inicio
    geom_point(data = sim_results[sum(!is.na(sim_results$NA_pop)), ], color = "darkgreen", size = 2) + # Punto final
    labs(title = paste("Gráfico de Fase Poblacional:", interaction_label),
         x = "Tamaño Poblacional de la Especie A",
         y = "Tamaño Poblacional de la Especie B") +
    theme_minimal()
  
  print(phase_plot)
} else {
  message("No se pudo generar el gráfico de fase ya que la simulación terminó prematuramente.")
}

# Inspirado en: León. J.A. 1974. Selection in context of interspecific competition. The American Naturalist, 108: 739-757.
