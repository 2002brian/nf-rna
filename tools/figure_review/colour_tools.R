# Colour metrics shared by the palette and scale review scripts (review aids only): CIELAB, CIEDE2000, colour-vision-
# deficiency simulation (colorspace) and WCAG contrast against white.  Requires the colorspace package.
lab <- function(hex) { m <- coords(as(hex2RGB(hex), "LAB")); matrix(m, ncol = 3, dimnames = list(hex, c("L", "a", "b"))) }
# CIEDE2000 (Sharma, Wu & Dalal 2005), vectorised over rows of two Lab matrices.
de2000 <- function(l1, l2) {
  rad <- pi / 180; deg <- 180 / pi
  C1 <- sqrt(l1[, 2]^2 + l1[, 3]^2); C2 <- sqrt(l2[, 2]^2 + l2[, 3]^2); Cb <- (C1 + C2) / 2
  G <- 0.5 * (1 - sqrt(Cb^7 / (Cb^7 + 25^7)))
  a1 <- (1 + G) * l1[, 2]; a2 <- (1 + G) * l2[, 2]
  C1p <- sqrt(a1^2 + l1[, 3]^2); C2p <- sqrt(a2^2 + l2[, 3]^2)
  h1 <- (atan2(l1[, 3], a1) * deg) %% 360; h2 <- (atan2(l2[, 3], a2) * deg) %% 360
  dL <- l2[, 1] - l1[, 1]; dC <- C2p - C1p
  dh <- h2 - h1; dh <- ifelse(C1p * C2p == 0, 0, ifelse(dh > 180, dh - 360, ifelse(dh < -180, dh + 360, dh)))
  dH <- 2 * sqrt(C1p * C2p) * sin(dh * rad / 2)
  Lb <- (l1[, 1] + l2[, 1]) / 2; Cpb <- (C1p + C2p) / 2
  hs <- h1 + h2; hb <- ifelse(C1p * C2p == 0, hs, ifelse(abs(h1 - h2) <= 180, hs / 2, ifelse(hs < 360, (hs + 360) / 2, (hs - 360) / 2)))
  T <- 1 - 0.17 * cos((hb - 30) * rad) + 0.24 * cos(2 * hb * rad) + 0.32 * cos((3 * hb + 6) * rad) - 0.20 * cos((4 * hb - 63) * rad)
  dTheta <- 30 * exp(-((hb - 275) / 25)^2); Rc <- 2 * sqrt(Cpb^7 / (Cpb^7 + 25^7))
  Sl <- 1 + 0.015 * (Lb - 50)^2 / sqrt(20 + (Lb - 50)^2); Sc <- 1 + 0.045 * Cpb; Sh <- 1 + 0.015 * Cpb * T
  Rt <- -sin(2 * dTheta * rad) * Rc
  sqrt((dL / Sl)^2 + (dC / Sc)^2 + (dH / Sh)^2 + Rt * (dC / Sc) * (dH / Sh))
}
pair_table <- function(hex) {
  idx <- which(upper.tri(matrix(0, length(hex), length(hex))), arr.ind = TRUE)
  L <- lab(hex); data.frame(i = idx[, 1], j = idx[, 2], de = de2000(L[idx[, 1], , drop = FALSE], L[idx[, 2], , drop = FALSE]))
}
visions <- list(normal = identity, deuteranopia = function(h) colorspace::deutan(h, severity = 1),
                protanopia = function(h) colorspace::protan(h, severity = 1), tritanopia = function(h) colorspace::tritan(h, severity = 1))
contrast_white <- function(hex) {
  lin <- function(v) ifelse(v <= 0.04045, v / 12.92, ((v + 0.055) / 1.055)^2.4)
  rgb <- t(col2rgb(hex)) / 255; y <- 0.2126 * lin(rgb[, 1]) + 0.7152 * lin(rgb[, 2]) + 0.0722 * lin(rgb[, 3])
  1.05 / (y + 0.05)
}
min_pair <- function(hex, vision) { p <- pair_table(visions[[vision]](hex)); k <- which.min(p$de); list(de = p$de[[k]], pair = paste(hex[p$i[[k]]], hex[p$j[[k]]])) }
cross_min <- function(a, b, vision) {
  La <- lab(visions[[vision]](a)); Lb <- lab(visions[[vision]](b))
  grid <- expand.grid(i = seq_along(a), j = seq_along(b)); de <- de2000(La[grid$i, , drop = FALSE], Lb[grid$j, , drop = FALSE]); k <- which.min(de)
  list(de = de[[k]], pair = paste(a[grid$i[[k]]], b[grid$j[[k]]]))
}
