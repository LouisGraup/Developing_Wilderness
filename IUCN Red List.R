
library(rredlist)
library(tidyverse)
library(data.table)

key='22e39dd150d9c1e4651b3e9134bc606169e0cb50f1958a9278a76f9582dc2dc8'

species = rl_sp(page=0, key=key)

df_species = species$result

for (i in 1:9) {

  out = rl_sp(page=i, key=key)

  df_species = rbind(df_species, out$result)

}

write.csv(df_species, file="species_list.csv")

regions = rl_regions(key)
regions = regions$results$identifier


for (i in 1:dim(df_species)[1]) {
  
  for (j in 1:length(regions)) {
  
    history = rl_history(name=df_species$scientific_name[i], key=key, region=regions[j])
    if (length(history$result)) {
      
      tempdf = history$result
      tempdf$Species = history$name
      tempdf$Class = df_species$class_name[i]
      tempdf$Region = history$region_identifier
      
      if (i==1) {
        dfhistory = tempdf
      }
      else {
        dfhistory = rbind(dfhistory, tempdf)
      }
      
    }
    
  }
  
}


for (i in 1:length(dfhistory$Species)) {
  
  dfhistory$count[i] = sum(dfhistory$Species[i]==dfhistory$Species)
  dfhistory$code_mod[i] = iucn_code(dfhistory$code[i])
  dfhistory$risk[i] = iucn_risk(dfhistory$code[i])
  
}

write.csv(dfhistory, file="iucn_historical_raw.csv")

dfManyHistories = filter(dfhistory, count>1 & code_mod!="NA")
all_years = sort(as.numeric(unique(dfManyHistories$year)))
all_species = unique(dfManyHistories$Species)

ggplot(dfManyHistories, aes(x=year, y=..count..)) + geom_bar(aes(fill=code_mod))

countries = rl_countries(key)
df_countries = countries$results

num_countries = dim(df_countries)[1]
dfSpeciesbyCountry = array(list(), num_countries)

for (i in 1:num_countries) {
  
  SpeciesbyCountry = rl_sp_country(df_countries$isocode[i], key=key)
  
  # loop through each country's species and grab all years of record
  
  if (SpeciesbyCountry$count>0) {
  
  dfSpeciesbyCountry[[i]] = list(country=df_countries$isocode[i], history=filter(dfManyHistories, Species==SpeciesbyCountry$result$scientific_name[1]))
  
    for (j in 2:SpeciesbyCountry$count) {
      dfSpeciesbyCountry[[i]]$history = rbind(dfSpeciesbyCountry[[i]]$history, filter(dfManyHistories, Species==SpeciesbyCountry$result$scientific_name[j]))
    }
  }
  else {
    dfSpeciesbyCountry[[i]]$country = df_countries$isocode[i]
  }
}

save(dfSpeciesbyCountry, file="SpeciesbyCountry.Rda")

dfCountryRisk = data.frame(matrix(nrow=num_countries, ncol=length(all_years)+1))
colnames(dfCountryRisk) = c("Country",all_years)

for (i in 1:num_countries) {
  
  CountryRisk = dfSpeciesbyCountry[[i]]$history
  CountryRisk$year = as.numeric(CountryRisk$year)
  dfCountryRisk$Country[i] = dfSpeciesbyCountry[[i]]$country
  
  if (length(CountryRisk$year)>0) {
  
  years = unique(CountryRisk$year)
  
  for (j in 1:length(years)) {
    # loop through all years of record and calculate RLI
    CountryRiskbyYear = filter(CountryRisk, year==years[j])
    dfCountryRisk[[which(all_years==years[j])+1]][i] = iucn_rli(CountryRiskbyYear$risk)
    
  }
  
  }
  
}

save(dfCountryRisk, file="CountryRiskbyYear.Rda")

# get list of countries for each species

for (i in 1:length(all_species)) {
  
  CountrybySpecies = rl_occ_country(name=all_species[i], key=key)
  CountrybySpecies$result$Species = all_species[i]
  dtCountrybySpecies = setDT(CountrybySpecies$result)
  
  # combine with histories to obtain risk level for each country, in each year
  tempHistory = setDT(filter(dfManyHistories, Species==all_species[i]))
  dtHistoriesbyCountry = tempHistory[dtCountrybySpecies, on="Species", allow.cartesian=TRUE]
  
  if (i==1) {
    dfCountrybySpecies = dtHistoriesbyCountry
  }
  else {
    dfCountrybySpecies = rbind(dfCountrybySpecies, dtHistoriesbyCountry, fill=TRUE)
  }
  
}

save(dfCountrybySpecies, file="CountrybySpecies.Rda")

# unique species ID list
species_id = array(0, length(all_species))

for (i in 1:length(all_species)) {
  
  species_id[i] = df_species$taxonid[which(df_species$scientific_name==all_species[i])]
  
}


iucn_rli = function(risk) {
  
  # calculate Red List Index as defined by IUCN
  RLI = 1 - (sum(risk) / (length(risk)*5))
  return(RLI)
  
}



iucn_risk = function(cd) {
  
  if (is.na(cd)) {
    return(0)
  }
  
  
  if (cd == "LC" || cd == "LR/lc" || cd == "LR/cd" || cd == "LR/nt" || cd == "O") {
    risk = 0
  }
  else if (cd == "NT" || cd == "nt" || cd == "R") {
    risk = 1
  }
  else if (cd == "VU" || cd == "V") {
    risk = 2
  }
  else if (cd == "E" || cd == "EN" || cd == "T") {
    risk = 3
  }
  else if (cd == "CR") {
    risk = 4
  }
  else if (cd == "EX" || cd == "Ex" || cd == "Ex?" || cd == "EW" || cd == "RE") {
    risk = 5
  }
  else {
    risk = 0
  }
  
  return(risk)
  
}

iucn_code = function(cd) {
  
  if (is.na(cd)) {
    return("NA")
  }
  
  if (cd == "LC" || cd == "LR/lc" || cd == "LR/cd" || cd == "LR/nt" || cd == "O") {
    code = "LC"
  }
  else if (cd == "NT" || cd == "nt" || cd == "R") {
    code = "NT"
  }
  else if (cd == "VU" || cd == "V") {
    code = "VU"
  }
  else if (cd == "E" || cd == "EN" || cd == "T") {
    code = "EN"
  }
  else if (cd == "CR") {
    code = "CR"
  }
  else if (cd == "EX" || cd == "Ex" || cd == "Ex?" || cd == "EW") {
    code = "EX"
  }
  else if (cd == "RE") {
    code = "RE"
  }
  else {
    code = "NA"
  }
  
  return(code)
  
}
