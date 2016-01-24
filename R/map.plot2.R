map.plot2 <- function (data, trait = NULL, trait.scale = "same", col.chr = NULL, 
                       col.trait = NULL, type = "hist", cex = 0.4, lwd = 1, cex.axis = 0.4, 
                       cex.trait = 0.8, jump = 5) 
{
  transp <- function(col, alpha = 0.5) {
    res <- apply(col2rgb(col), 2, function(c) rgb(c[1]/255, 
                                                  c[2]/255, c[3]/255, alpha))
    return(res)
  }
  len <- numeric()
  for (i in 1:max(unique(data$LG))) {
    len[i] <- max(data[which(data$LG == i), "Position"])
  }
  coree <- which(len == max(len))
  coree1 <- data[which(data$LG == coree), ]
  linesss <- coree1$Position/max(coree1$Position)
  if (!is.null(trait)) {
    cols <- 1
  }
  else {
    cols <- 0
  }
  extra <- length(len) + cols * length(len)
  fact <- 1/extra
  fact2 <- fact + (fact * cols)
  if (!is.null(col.trait)) {
    col.trait <- col.trait
  }
  else {
    col.trait <- c(1:6, 1:6)
  }
  plot.new()
  for (j in 1:max(unique(data$LG))) {
    prov <- data[which(data$LG == j), ]
    dddd <- prov[1, ]
    dddd2 <- prov[dim(prov)[1], ]
    dddd[1, which(names(dddd) != "LG")] <- 0
    dddd2[1, which(names(dddd) != "LG" & names(dddd) != "Position")] <- 0
    prov <- rbind(dddd, prov, dddd2)
    chr <- prov$Position/max(coree1$Position)
    ruler <- 1 - (c(seq(0, max(prov$Position), by = jump), 
                    round(max(prov$Position), 0))/max(coree1$Position))
    ruler2 <- c(seq(0, max(prov$Position), by = jump), round(max(prov$Position), 
                                                             0))
    if (!is.null(trait)) {
      sss <- (fact2 * j) - fact
    }
    else {
      sss <- (fact2 * j)
    }
    dd2 <- density(chr, n = length(chr))$y
    dd <- sort(density(chr, n = length(chr))$y, decreasing = T)
    if (!is.null(col.chr)) {
      hc <- colorRampPalette(c(col.chr[1], col.chr[2]))(length(dd))
    }
    else {
      hc <- gray.colors(n = length(dd), start = 0, end = 0.6, 
                        gamma = 2.2, alpha = NULL)
    }
    for (k in 1:length(dd2)) {
      ooo <- which(dd == dd2[k])
      lines(y = c(1 - chr[k], 1 - chr[k]), x = c(sss, sss - 
                                                   (fact/3)), lwd = lwd, col = hc[ooo])
    }
    lines(y = c(1, 1 - max(chr)), x = c(sss, sss), lwd = 3)
    lines(y = c(1, 1 - max(chr)), x = c(sss - (fact/3), sss - 
                                          (fact/3)), lwd = 3)
    text(x = sss - (fact/1.6), y = ruler, labels = ruler2, 
         cex = cex)
    axis(3, at = (sss - (fact/3)), labels = paste("LG", j, 
                                                  sep = ""), cex.axis = cex.axis, font = 2)
    plotrix::draw.arc((sss + sss - (fact/3))/2, 1 - max(chr), 
                      (sss - (sss - (fact/3)))/2, deg1 = 180, deg2 = 360, 
                      col = "black", lwd = 2, lend = 1)
    plotrix::draw.arc((sss + sss - (fact/3))/2, 1, (sss - 
                                                      (sss - (fact/3)))/2, deg1 = 0, deg2 = 180, col = "black", 
                      lwd = 2, lend = 1)
    if (!is.null(trait)) {
      if (is.numeric(prov[, trait])) {
        w1 <- which(names(data) == trait)
        if (trait.scale == "same") {
          bobo <- max(data[, trait], na.rm = TRUE)
        }
        else {
          bobo <- max(prov[, trait], na.rm = TRUE)
        }
        dotss <- fact * (prov[, trait]/bobo)
        dotss2 <- sss + (fact/8) + dotss
        sections <- (bobo - min(prov[, trait], na.rm = TRUE))/5
        sections2 <- seq(min(prov[, trait], na.rm = TRUE), 
                         bobo, by = sections)
        sections3 <- fact * (sections2/bobo)
        sections4 <- sss + (fact/8) + sections3
        for (d in 1:length(sections4)) {
          lines(x = c(sections4[d], sections4[d]), y = c(1, 
                                                         1 - max(chr)), col = "black", lty = 3, lwd = 0.5)
          text(x = sections4[d], y = 1, labels = round(sections2[d], 
                                                       1), cex = 0.4, srt = 270)
        }
        if (type == "dot") {
          points(y = 1 - chr, x = dotss2, pch = 20, cex = cex.trait, 
                 col = transp(col.trait[j], 0.6))
        }
        if (type == "line") {
          polygon(y = 1 - chr, x = dotss2, pch = 20, 
                  cex = cex.trait, col = transp(col.trait[j], 
                                                0.4))
          lines(y = 1 - chr, x = dotss2, pch = 20, cex = cex.trait, 
                col = transp(col.trait[j], 0.6))
        }
        if (type == "hist") {
          for (l in 1:length(dotss2)) {
            lines(x = c(sss + (fact/8), dotss2[l]), y = c(1 - 
                                                            chr[l], 1 - chr[l]), lwd = cex.trait, col = transp(col.trait[j], 
                                                                                                               0.8))
          }
        }
        axis(3, at = sss + (fact/2), labels = trait, 
             cex.axis = cex.axis)
      }
      if (is.factor(prov[, trait])) {
        riel <- sss + (fact/2)
        lines(x = c(riel, riel), y = c(1, 1 - max(chr)), 
              col = transp(col.trait[j], 0.8))
        ww2 <- which(!is.na(prov[, trait]))
        if (length(ww2) > 0) {
          yy <- 1 - chr[ww2]
          points(y = yy, x = rep(riel, length(yy)), pch = 15, 
                 cex = 0.9, col = transp(col.trait[j], 0.6))
          text(x = riel + (fact/2), y = yy, labels = prov[ww2, 
                                                          trait], cex = 0.3)
        }
        axis(3, at = riel, labels = trait, cex.axis = cex)
      }
      if (is.character(prov[, trait])) {
        riel <- sss + (fact/1.8)
        ww2 <- which(!is.na(prov[, trait]))
        if (length(ww2) > 0) {
          yy <- 1 - chr[ww2]
          text(x = riel, y = yy, labels = prov[ww2, trait], 
               cex = 0.17)
        }
        axis(3, at = riel, labels = trait, cex.axis = cex.axis)
      }
    }
  }
}