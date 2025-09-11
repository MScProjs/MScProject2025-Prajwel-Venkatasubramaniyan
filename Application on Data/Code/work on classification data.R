preprocess <- function(data)
{
  X <- data[,-c(1,2,4)]
  y <- data[,4]

  dat <- as.data.frame(X)
  dat$target <- y

  return(dat)
}


alphavirus <- preprocess(read.csv("Alphavirus.csv"))
instance <- list(alias = "Alphavirus", Data = alphavirus)
algorithms <- c("svm", "random_forest")
r9 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r9)
plot.nreps(r9)

bordetella <- preprocess(read.csv("Bordetella pertussis.csv"))
instance <- list(alias = "Bordetella pertussis", Data = bordetella)
algorithms <- c("svm", "random_forest")
r8 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r8)
plot.nreps(r8)

campylobactor <- preprocess(read.csv("Campylobacter jejuni.csv"))
instance <- list(alias = "Campylobacter jejuni", Data = campylobactor)

chlamydia <- preprocess(read.csv("Chlamydia trachomatis.csv"))
instance <- list(alias = "Chlamydia trachomatis", Data = chlamydia)

chostridioides <- preprocess(read.csv("Clostridioides difficile.csv"))
instance <- list(alias = "Clostridioides difficile", Data = chostridioides)

dengue <- preprocess(read.csv("Dengue virus.csv"))
instance <- list(alias = "Dengue virus", Data = dengue)

enterovirus <- preprocess(read.csv("Enterovirus C.csv"))
instance <- list(alias = "Enterovirus C", Data = enterovirus)
algorithms <- c("svm", "random_forest")
r10 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r10)
plot.nreps(r10)

escherichia <- preprocess(read.csv("Escherichia coli.csv"))
instance <- list(alias = "Escherichia coli", Data = escherichia)

haemophilus <- preprocess(read.csv("Haemophilus influenzae.csv"))
instance <- list(alias = "Haemophilus influenzae", Data = haemophilus)

hepacivirus <- preprocess(read.csv("Hepacivirus hominis.csv"))
instance <- list(alias = "Hepacivirus hominis", Data = hepacivirus)

human <- preprocess(read.csv("human gammaherpesvirus 4.csv"))
instance <- list(alias = "human gammaherpesvirus 4", Data = human)


Human <- preprocess(read.csv("Human immunodeficiency virus 1.csv"))
instance <- list(alias = "Human immunodeficiency virus ", Data = Human)
algorithms <- c("svm", "random_forest")
r1 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r1)
plot.nreps(r1)

influenza <- preprocess(read.csv("Influenza A virus.csv"))
instance <- list(alias = "influenza", Data = influenza)

klebsiella <- preprocess(read.csv("Klebsiella pneumoniae.csv"))
instance <- list(alias = "klebsiella pneumoniae", Data = klebsiella)
algorithms <- c("svm", "random_forest")
r2 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r2)
plot.nreps(r2)

leishmania <- preprocess(read.csv("Leishmania infantum.csv"))
instance <- list(alias = "Leishmania infantum", Data = leishmania)

leptospira <- preprocess(read.csv("Leptospira.csv"))
instance <- list(alias = "Leptospira", Data = leptospira)
algorithms <- c("svm", "random_forest")
r3 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r3)
plot.nreps(r3)

measles <- preprocess(read.csv("Measles morbillivirus.csv"))
instance <- list(alias = "Measles morbillivirus", Data = measles)

mycobacterium <- preprocess(read.csv("Mycobacterium tuberculosis.csv"))
instance <- list(alias = "Mycobacterium tuberculosis", Data = mycobacterium)

neisseria <- preprocess(read.csv("Neisseria gonorrhoeae.csv"))
instance <- list(alias = "Neisseria gonorrhoeae", Data = neisseria)
algorithms <- c("svm", "random_forest")
r4 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r4)
plot.nreps(r4)

onchocerca <- preprocess(read.csv("Onchocerca volvulus.csv"))
instance <- list(alias = "Onchocerca volvulus", Data = onchocerca)

orthoebolavirus <- preprocess(read.csv("Orthoebolavirus zairense.csv"))
instance <- list(alias = "Orthoebolavirus zairense", Data = orthoebolavirus)

orthopoxvirus <- preprocess(read.csv("Orthopoxvirus.csv"))
instance <- list(alias = "Orthopoxvirus", Data = orthopoxvirus)
algorithms <- c("svm", "random_forest")
r5 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r5)
plot.nreps(r5)

plasmodium <- preprocess(read.csv("Plasmodium falciparum.csv"))
instance <- list(alias = "Plasmodium falciparum", Data = plasmodium)

pseudomonas <- preprocess(read.csv("Pseudomonas aeruginosa.csv"))
instance <- list(alias = "Pseudomonas aeruginosa", Data = pseudomonas)
algorithms <- c("svm", "random_forest")
r6 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r6)
plot.nreps(r6)

schistosoma <- preprocess(read.csv("Schistosoma mansoni.csv"))
instance <- list(alias = "Schistosoma mansoni", Data = schistosoma)

severe <- preprocess(read.csv("Severe acute respiratory syndrome-related coronavirus.csv"))
instance <- list(alias = "Severe", Data = severe)

staphylococcus <- preprocess(read.csv("Staphylococcus aureus.csv"))
instance <- list(alias = "Staphylococcus aureus", Data = staphylococcus)

Streptococcus <- preprocess(read.csv("Streptococcus pyogenes.csv"))
instance <- list(alias = "Streptococcus pyogenes", Data = Streptococcus)


toxoplasma <- preprocess(read.csv("Toxoplasma gondii.csv"))
instance <- list(alias = "Toxoplasma gondii", Data = toxoplasma)
algorithms <- c("svm", "random_forest")
r7 <- calc_nreps_classification(instance, algorithms)
summary.nreps(r7)
plot.nreps(r7)

yellow <- preprocess(read.csv("Yellow fever virus.csv"))
instance <- list(alias = "Yellow fever virus", Data = yellow)

#####################################################################################

instances <- list(list(alias = "Alphavirus", Data = alphavirus),
                  list(alias = "Bordetella pertussis", Data = bordetella),
                  list(alias = "Campylobacter jejuni", Data = campylobactor),
                  list(alias = "Chlamydia trachomatis", Data = chlamydia),
                  list(alias = "Clostridioides difficile", Data = chostridioides),

                  list(alias = "Dengue virus", Data = dengue),
                  list(alias = "Enterovirus C", Data = enterovirus),
                  list(alias = "Escherichia coli", Data = escherichia),
                  list(alias = "Haemophilus influenzae", Data = haemophilus),
                  list(alias = "Hepacivirus hominis", Data = hepacivirus),

                  list(alias = "human gammaherpesvirus 4", Data = human),
                  list(alias = "Human immunodeficiency virus ", Data = Human),
                  list(alias = "influenza", Data = influenza),
                  list(alias = "klebsiella pneumoniae", Data = klebsiella),
                  list(alias = "Leishmania infantum", Data = leishmania),

                  list(alias = "Leptospira", Data = leptospira),
                  list(alias = "Measles morbillivirus", Data = measles),
                  list(alias = "Mycobacterium tuberculosis", Data = mycobacterium),
                  list(alias = "Neisseria gonorrhoeae", Data = neisseria),
                  list(alias = "Onchocerca volvulus", Data = onchocerca),

                  list(alias = "Orthoebolavirus zairense", Data = orthoebolavirus),
                  list(alias = "Orthopoxvirus", Data = orthopoxvirus),
                  list(alias = "Plasmodium falciparum", Data = plasmodium),
                  list(alias = "Pseudomonas aeruginosa", Data = pseudomonas),
                  list(alias = "Schistosoma mansoni", Data = schistosoma),

                  list(alias = "Severe", Data = severe),
                  list(alias = "Staphylococcus aureus", Data = staphylococcus),
                  list(alias = "Streptococcus pyogenes", Data = Streptococcus),
                  list(alias = "Toxoplasma gondii", Data = toxoplasma),
                  list(alias = "Yellow fever virus", Data = yellow)

)
algorithms <- c("svm", "random_forest")
final_result_classification <- run_experiment_classification(instances, algorithms)
