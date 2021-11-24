
library(conflicted) # A very convenient package to avoid and resolve namespace conflicts
library(janitor) # A package to clean tables
library(tidyr) # A reference package to manipulate data
library(dplyr) # A reference package to manipulate data
library(ggplot2) # A reference package to plot data
library(glue) # A package to put variables in strings.

conflict_prefer("filter", "dplyr") # Conflict resolution for the 'filter' function 

theme_set(theme_bw(base_size = 13)) # Theme to use for all plots

qtlDF <- read.table(
    file = "C:\\Projects\\IBIP21\\idx3005isCausal.csv.gz",
    header = T,
    stringsAsFactors = F,
    sep = ","
) %>% 
    clean_names() %>% 
    rename(
        chromosome = chr,
        position = pos,
        genotype = gt,
        phenotype = pheno,
        sample = subid
    ) %>% 
    select(
        -x
    )

head(qtlDF)

snpDF <- qtlDF %>% 
    filter(
        position == 37000661
    )

table(
    snpDF$genotype
)

snpLmResults <- lm(
    formula = phenotype ~ genotype,
    data = snpDF
)

summary(snpLmResults)

ggplot(
    data = snpDF
    ) +
    geom_point(
        mapping = aes(
            x = genotype,
            y = phenotype
        ),
        alpha = 0.05
    ) + 
    geom_smooth(
        mapping = aes(
            x = genotype,
            y = phenotype
        ),
        method = "lm"
    )

qtl1DF <- qtlDF %>% 
    filter(
        position >= 37500000 & position <= 38100000
    ) %>% 
    mutate(
        protein_A = rnorm(n = n()),
        protein_B = phenotype,
        protein_C = rnorm(n = n()),
        protein_D = rnorm(n = n())
    )

write.table(
    x = qtl1DF,
    file = gzfile("pages/proteogenomics/pQTL/resources/pqtl.gz"),
    col.names = T,
    row.names = F,
    sep = "\t",
    quote = F
)

lmResults <- data.frame(
    position = unique(qtl1DF$pos),
    stringsAsFactors = F
)

lmResults$p_A <- sapply(
    X = lmResults$position, 
    FUN = function (position) summary(
        lm(
            formula = protein_A ~ genotype,
            data = qtl1DF[qtl1DF$position == position,]
        )
    )$coef["genotype", "Pr(>|t|)"]
)

lmResults$p_B <- sapply(
    X = lmResults$position, 
    FUN = function (position) summary(
        lm(
            formula = protein_B ~ genotype,
            data = qtl1DF[qtl1DF$position == position,]
        )
    )$coef["genotype", "Pr(>|t|)"]
)

lmResults$p_C <- sapply(
    X = lmResults$position, 
    FUN = function (position) summary(
        lm(
            formula = protein_C ~ genotype,
            data = qtl1DF[qtl1DF$position == position,]
        )
    )$coef["genotype", "Pr(>|t|)"]
)

lmResults$p_D <- sapply(
    X = lmResults$position, 
    FUN = function (position) summary(
        lm(
            formula = protein_D ~ genotype,
            data = qtl1DF[qtl1DF$position == position,]
        )
    )$coef["genotype", "Pr(>|t|)"]
)

plot(lmResults$position, -log10(lmResults$p_A))
plot(lmResults$position, -log10(lmResults$p_B))
plot(lmResults$position, -log10(lmResults$p_C))
plot(lmResults$position, -log10(lmResults$p_D))

plotDF <- lmResults %>% 
    pivot_longer(
        cols = c("p_A", "p_B", "p_C", "p_D"),
        names_prefix = "p_"
    )

ggplot(
    data = plotDF
) +
    geom_hline(
        yintercept = 0
    ) +
    geom_point(
        mapping = aes(
            x = position,
            y = -log10(value),
            col = name
        )
    ) +
    scale_y_continuous(
        name = "p-value [-log10]"
    ) +
    facet_grid(
        name ~ .
    ) +
    theme(
        legend.position = "none"
    )


getCiUp <- function(n, confidence) {
    
    return(
        qbeta(p = (1-confidence)/2, shape1 = 1:n, shape2 = n - 1:n + 1)
    )
}

getCiDown <- function(n, confidence) {
    
    return(
        qbeta(p = 1-(1-confidence)/2, shape1 = 1:n, shape2 = n - 1:n + 1)
    )
}

confidence <- 0.95

qqPlotDF <- plotDF %>% 
    group_by(
        name
    ) %>% 
    arrange(
        value
    ) %>% 
    mutate(
        logP = -log10(value),
        expectedLogP = -log10(sort(ppoints(n = n()))),
        ciUp = getCiUp(
            n = n(), 
            confidence = confidence
        ),
        ciDown = getCiDown(
            n = n(), 
            confidence = confidence
        ),
        ciUpLog = -log10(ciDown),
        ciDownLog = -log10(ciUp)
    )

ggplot(
    data = qqPlotDF
) +
    geom_ribbon(
        data = qqPlotDF,
        mapping = aes(
            x = expectedLogP,
            ymin = ciDownLog,
            ymax = ciUpLog
        ),
        fill = "black",
        alpha = 0.2
    ) +
    geom_abline(
        linetype = "dotted"
    ) +
    geom_point(
        data = test,
        mapping = aes(
            x = expectedLogP,
            y = logP,
            col = name
        )
    ) + 
    scale_x_continuous(
        name = "Expected p-value [-log10]"
    ) + 
    scale_y_continuous(
        name = "Observed p-value [-log10]"
    ) +
    theme(
        legend.position = "top"
    )

