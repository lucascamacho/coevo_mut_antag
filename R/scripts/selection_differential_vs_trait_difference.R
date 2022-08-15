# Plotting selection differentials

# cleaning wd and loading packages ------------------------------ 
rm(list = ls())
set.seed(10)
if(!require(ggplot2)) {install.packages("ggplot2"); library(ggplot2)}
if(!require(latex2exp)) {install.packages("latex2exp"); library(latex2exp)}
if(!require(gridExtra)) {install.packages("gridExtra"); library(gridExtra)}

# environmental selection differential ------------------------------ 
theta <- 10
z_i <- rnorm(1000, 10, 2)
z_j <- rnorm(1000, 10, 2)
S_e <- theta - z_i
df <- data.frame(z_i = z_i, S_e = S_e)
p_e = ggplot(data = df, aes(x = z_i, y = S_e)) +
  geom_line(size = 1.5) +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab(TeX("Trait ($Z_i$)")) +
  ylab(TeX("Environmental selection differential ($S_e$)")) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# mutualistic selection differential ------------------------------ 
alpha <- 0.1
z_diff <- z_j - z_i
q <- exp(-alpha * (z_diff)^2)
S_m <- z_diff * q
df <- data.frame(z_diff = z_diff, S_m = S_m)
p_m = ggplot(data = df, aes(x = z_diff, y = S_m)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab(TeX("Trait difference ($Z_j - Z_i$)")) +
  ylab(TeX("Mutualistic selection differential ($S_m$)")) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# antagonistic selection differential ------------------------------ 
epsilon <- 5
delta <- ifelse(abs(z_diff) < epsilon, 1, 0)
sign <- ifelse(z_diff < 0, 1, -1)
S_a <- delta * q * (z_j + sign * epsilon - z_i)
df <- data.frame(z_diff = z_diff, S_a = S_a)
p_c = ggplot(data = df, aes(x = z_diff, y = S_a)) +
  geom_line(size = 1.5) +
  geom_vline(xintercept = 0, linetype = "dashed") +
  geom_hline(yintercept = 0, linetype = "dashed") +
  xlab(TeX("Trait difference ($Z_j - Z_i$)")) +
  ylab(TeX("Antagonistic selection differential ($S_a$)")) +
  theme_bw() +
  theme(plot.background = element_rect(fill = "transparent",colour = NA),
        panel.grid.major = element_blank(), panel.grid.minor = element_blank(),
        panel.background = element_blank(), axis.line = element_line(colour = "black"),
        axis.text.x = element_text(size = 13),
        axis.text.y = element_text(size = 13),
        axis.title = element_text(size = 16), 
        legend.key.size = unit(0.9, "cm"),
        legend.text = element_text(size = 13))

# final plot --------------------------------------------------------
plot_final = grid.arrange(p_m, p_c, p_e, nrow = 1)

ggsave(plot_final, filename = "Figura2.png", dpi = 600,
       width = 27, height = 16, units = "cm",  bg = "transparent")
