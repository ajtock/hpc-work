#promoters() in GenomicRanges 
function (x, upstream = 2000, downstream = 200, ...) 
{
    if (!isSingleNumber(upstream)) 
        stop("'upstream' must be a single integer")
    if (!is.integer(upstream)) 
        upstream <- as.numeric(upstream)
    if (!isSingleNumber(downstream)) 
        stop("'downstream' must be a single integer")
    if (!is.integer(downstream)) 
        downstream <- as.numeric(downstream)
    if (upstream < 0 | downstream < 0) 
        stop("'upstream' and 'downstream' must be integers >= 0")
    if (any(strand(x) == "*")) 
        warning("'*' ranges were treated as '+'")
    on_plus <- which(strand(x) == "+" | strand(x) == "*")
    on_plus_TSS <- start(x)[on_plus]
    start(x)[on_plus] <- on_plus_TSS - upstream
    end(x)[on_plus] <- on_plus_TSS + downstream - 1L
    on_minus <- which(strand(x) == "-")
    on_minus_TSS <- end(x)[on_minus]
    end(x)[on_minus] <- on_minus_TSS + upstream
    start(x)[on_minus] <- on_minus_TSS - downstream + 1L
    x
}

##Adapted from promoters() in GenomicRanges to extract regions relative to TTS
TTSplus <- function (x, upstream = 100, downstream = 1000, ...)
{
    if (!isSingleNumber(upstream))
        stop("'upstream' must be a single integer")
    if (!is.integer(upstream))
        upstream <- as.numeric(upstream)
    if (!isSingleNumber(downstream))
        stop("'downstream' must be a single integer")
    if (!is.integer(downstream))
        downstream <- as.numeric(downstream)
    if (downstream < 0)
        stop("'downstream' must be an integer >= 0")
#    if (upstream < 0 | downstream < 0)
#        stop("'upstream' and 'downstream' must be integers >= 0")
    if (any(strand(x) == "*"))
        warning("'*' ranges were treated as '+'")
    on_plus <- which(strand(x) == "+" | strand(x) == "*")
    on_plus_TTS <- end(x)[on_plus]
    start(x)[on_plus] <- on_plus_TTS - upstream
    end(x)[on_plus] <- on_plus_TTS + downstream
    on_minus <- which(strand(x) == "-")
    on_minus_TTS <- start(x)[on_minus]
    end(x)[on_minus] <- on_minus_TTS + upstream
    start(x)[on_minus] <- on_minus_TTS - downstream
    x
}

