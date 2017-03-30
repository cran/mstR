testListMST<-function (list, type = "start") 
{
    argNames = ls()
    if (!is.list(list)) 
        res <- list(test = FALSE, message = paste(deparse(substitute(list)), 
            " is not a list", sep = ""))
    else {
        if (is.null(names(list))) 
            res <- list(test = FALSE, message = paste("list '", 
                deparse(substitute(list)), "' has no argument names", 
                sep = ""))
        else {
            elements <- switch(type, start = c("fixModule", "seed", 
                "theta", "D"), test = c("method", "priorDist", "priorPar", "range", 
                  "D", "parInt", "moduleSelect", "constantPatt", "cutoff"), 
                 final = c("method", "priorDist", "priorPar", "range", "D", "alpha", 
                  "parInt"))
            if (is.null(elements)) 
                res <- list(test = FALSE, message = paste("invalid 'type' argument ('", 
                  type, "' is not allowed)", sep = ""))
            else {
                if (length(list) > length(elements)) 
                  res <- list(test = FALSE, mesage = paste("too many elements in ", 
                    deparse(substitute(list)), " for type '", 
                    type, "'", sep = ""))
                else {
                  res <- list(test = TRUE, message = "ok")
                  i <- 0
                  repeat {
                    i <- i + 1
                    if (sum(names(list)[i] == elements) == 0) {
                      res$test <- FALSE
                      break
                    }
                    else {
                      if (i == length(list)) 
                        break
                    }
                  }
                  if (!res$test) {
                    texte <- switch(i, `1` = "st", `2` = "nd", 
                      `3` = "rd")
                    if (is.null(texte)) 
                      texte <- "th"
                    res$message <- paste("invalid name '", names(list)[i], 
                      "' for ", i, texte, " element of '", deparse(substitute(list)), 
                      "'", sep = "")
                  }
                  else {
                    intNames <- c("fixModule")
                    seedNames <- c("seed")
                    numNames <- c("alpha", "D", "theta")
                    metNames <- c("method")
                    priorNames <- c("priorDist")
                    parNames <- c("priorPar", "range")
                    eapNames <- c("parInt")
                    moduleSelectNames <- c("moduleSelect")
                    constantNames <- c("constantPatt")
                    matrixNames <- c("cutoff")
                    logicNames <- c("cb.control")
                    i <- 0
                    repeat {
                      i <- i + 1
                      vect <- c(sum(names(list)[i] == intNames), 
                        sum(names(list)[i] == seedNames), sum(names(list)[i] == 
                          numNames), sum(names(list)[i] == metNames), 
                          sum(names(list)[i] == priorNames), sum(names(list)[i] == 
                          parNames), sum(names(list)[i] == eapNames), 
                          sum(names(list)[i] == moduleSelectNames), 
                          sum(names(list)[i] == constantNames), 
                          sum(names(list)[i] == matrixNames), 
                          sum(names(list)[i] == logicNames))
                      ind <- (1:11)[vect == 1]
                      prov <- switch(ind, `1` = ifelse(is.null(list[[i]]), 
                        TRUE, ifelse(is.numeric(list[[i]]), ifelse(max(abs(list[[i]] - 
                          round(list[[i]]))) <= 1e-04, TRUE, 
                          FALSE), FALSE)), `2` = ifelse(is.null(list[[i]]), 
                        TRUE, ifelse(is.numeric(list[[i]]), ifelse(length(list[[i]]) == 
                          1, TRUE, FALSE), ifelse(is.na(list[[i]]), 
                          TRUE, FALSE))), `3` = (is.numeric(list[[i]]) & 
                        length(list[[i]]) == 1), `4` = (is.list(list[[i]]) == 
                        FALSE & length(list[[i]]) == 1 & sum(list[[i]] == 
                        c("ML", "BM", "WL", "EAP","score")) == 1), `5` = (is.list(list[[i]]) == 
                        FALSE & length(list[[i]]) == 1 & sum(list[[i]] == 
                        c("norm", "unif", "Jeffreys")) == 1), 
                        `6` = (is.numeric(list[[i]]) & length(list[[i]]) == 
                          2), `7` = (!is.list(list[[i]]) & 
                          length(list[[i]]) == 3 & is.numeric(list[[i]]) == 
                          TRUE & abs(list[[i]][3] - round(list[[i]][3])) <= 
                          1e-04), `8` = (is.list(list[[i]]) == 
                          FALSE & length(list[[i]]) == 1 & sum(list[[i]] == 
                          c("MFI", "MLMWI", "MPMWI", "MKL", "MKLP",  
                          "random")) == 1), `9` = ifelse(is.null(list[[i]]), 
                          TRUE, ifelse(sum(list[[i]] == c("fixed4", 
                          "fixed7", "var", "BM", "EAP", "WL")) == 
                          1, TRUE, FALSE)), `10` = ifelse(is.matrix(list[[i]]), 
                          TRUE, FALSE), `11` = ifelse(is.logical(list[[i]]), 
                          TRUE, FALSE))
                      if (!prov) {
                        res$test <- FALSE
                        res$message <- switch(ind, `1` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a vector of integer values or NULL", 
                          sep = ""), `2` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a single numeric value or NULL", 
                          sep = ""), `3` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a single numeric value", 
                          sep = ""), `4` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be either 'ML', 'BM', 'EAP', 'WL' or 'score'", 
                          sep = ""), `5` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be either 'norm', 'unif'or 'Jeffreys'", 
                          sep = ""), `6` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a vector of two numeric values", 
                          sep = ""), `7` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be a vector of two numeric and", 
                          "\n", " one integer components", sep = ""), 
                          `8` = paste("element '", names(list)[i], 
                          "' of '", deparse(substitute(list)), 
                          "' must be either 'MFI', 'MLMWI',", 
                          "\n", " 'MPMWI', 'MKL', 'MKLP' or 'random'", 
                          sep = ""), `9` = paste("element '", 
                          names(list)[i], "' of '", deparse(substitute(list)), 
                          "' must be either 'fixed4', 'fixed7', 'var' or NULL", 
                          sep = ""), `10` = paste("element '", 
                            names(list)[i], "' of '", deparse(substitute(list)), 
                            "' must be a suitable matrix of cut-off values", 
                            sep = ""), `11` = paste("element '", 
                            names(list)[i], "' of '", deparse(substitute(list)), 
                            "' must be either 'TRUE' or 'FALSE'", 
                            sep = ""))
                        break
                      }
                      else {
                        if (i == length(list)) 
                          break
                      }
                    }
                  }
                }
            }
        }
    }
    return(res)
}
