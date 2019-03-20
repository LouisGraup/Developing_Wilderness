
library(rredlist)
library(tidyverse)
library(data.table)

key='22e39dd150d9c1e4651b3e9134bc606169e0cb50f1958a9278a76f9582dc2dc8' #unique key for IUCN data access

# grab list of all species in IUCN database
species = rl_sp(page=0, key=key) # call to server

df_species = species$result

for (i in 1:9) {

  out = rl_sp(page=i, key=key) # call to server

  df_species = rbind(df_species, out$result)

}

write.csv(df_species, file="species_list.csv")

# loop through regions to grab historical classifications of each species

regions = rl_regions(key) # call to server
regions = regions$results$identifier


for (i in 1:dim(df_species)[1]) {
  
  for (j in 1:length(regions)) {
  
    history = rl_history(name=df_species$scientific_name[i], key=key, region=regions[j]) # call to server
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

# calculate # of classifications per species and IUCN Red List code/risk
for (i in 1:length(dfhistory$Species)) {
  
  dfhistory$count[i] = sum(dfhistory$Species[i]==dfhistory$Species)
  dfhistory$code_mod[i] = iucn_code(dfhistory$code[i])
  dfhistory$risk[i] = iucn_risk(dfhistory$code[i])
  
}

write.csv(dfhistory, file="iucn_historical_raw.csv")

# remove species with only one classification

dfManyHistories = filter(dfhistory, count>1 & code_mod!="NA")

# get list of years and species
all_years = sort(as.numeric(unique(dfManyHistories$year)))
all_species = unique(dfManyHistories$Species)

ggplot(dfManyHistories, aes(x=year, y=..count..)) + geom_bar(aes(fill=code_mod))

# get list of countries
countries = rl_countries(key)
df_countries = countries$results

num_countries = dim(df_countries)[1]

# old methodology
# dfSpeciesbyCountry = array(list(), num_countries)
# 
# for (i in 1:num_countries) {
#   
#   SpeciesbyCountry = rl_sp_country(df_countries$isocode[i], key=key)
#   
#   # loop through each country's species and grab all years of record
#   
#   if (SpeciesbyCountry$count>0) {
#   
#   dfSpeciesbyCountry[[i]] = list(country=df_countries$isocode[i], history=filter(dfManyHistories, Species==SpeciesbyCountry$result$scientific_name[1]))
#   
#     for (j in 2:SpeciesbyCountry$count) {
#       dfSpeciesbyCountry[[i]]$history = rbind(dfSpeciesbyCountry[[i]]$history, filter(dfManyHistories, Species==SpeciesbyCountry$result$scientific_name[j]))
#     }
#   }
#   else {
#     dfSpeciesbyCountry[[i]]$country = df_countries$isocode[i]
#   }
# }
# 
# save(dfSpeciesbyCountry, file="SpeciesbyCountry.Rda")

# new methodology
# get list of countries for each species

for (i in 1:length(all_species)) {
  
  CountrybySpecies = rl_occ_country(name=all_species[i], key=key) # call to server
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


gdp_species$weight = 1
# update data with gdp and habitat extent information

for (i in 1:length(gdp_species$binomial)) {
  # calculate weight for each species in each country 
  # based on percentage of species' habitat in an individual country over the total habitat extent
  gdp_species$weight[i] = gdp_species$area[i]/sum(filter(gdp_species, binomial==gdp_species$binomial[i])$area)
  
}

# identifier column
gdp_species$id = paste(gdp_species$country, gdp_species$binomial)

# combine dfCountrybySpecies and gdp_species based on identifier
for (i in 28994:length(dfCountrybySpecies$Species)) {
  
  j = match(paste(dfCountrybySpecies$country[i], dfCountrybySpecies$Species[i]), gdp_species$id)
  if (!is.na(j)) {
  dfCountrybySpecies$country_area_sqmil[i] = gdp_species$area_sqmil[j]
  dfCountrybySpecies$country_area_sqkm[i] = gdp_species$area_sqkm[j]
  dfCountrybySpecies$species_area[i] = gdp_species$area[j]
  dfCountrybySpecies$pct_species_country[i] = gdp_species$pct_specie_country[j]
  dfCountrybySpecies$yr2000[i] = gdp_species$yr2000[j]
  dfCountrybySpecies$yr2001[i] = gdp_species$yr2001[j]
  dfCountrybySpecies$yr2002[i] = gdp_species$yr2002[j]
  dfCountrybySpecies$yr2003[i] = gdp_species$yr2003[j]
  dfCountrybySpecies$yr2004[i] = gdp_species$yr2004[j]
  dfCountrybySpecies$yr2005[i] = gdp_species$yr2005[j]
  dfCountrybySpecies$yr2006[i] = gdp_species$yr2006[j]
  dfCountrybySpecies$yr2007[i] = gdp_species$yr2007[j]
  dfCountrybySpecies$yr2008[i] = gdp_species$yr2008[j]
  dfCountrybySpecies$yr2009[i] = gdp_species$yr2009[j]
  dfCountrybySpecies$yr2010[i] = gdp_species$yr2010[j]
  dfCountrybySpecies$yr2011[i] = gdp_species$yr2011[j]
  dfCountrybySpecies$yr2012[i] = gdp_species$yr2012[j]
  dfCountrybySpecies$yr2013[i] = gdp_species$yr2013[j]
  dfCountrybySpecies$yr2014[i] = gdp_species$yr2014[j]
  dfCountrybySpecies$yr2015[i] = gdp_species$yr2015[j]
  dfCountrybySpecies$weight[i] = gdp_species$weight[j]
  }
}

save(dfCountrybySpecies, file="CountrybySpecies.Rda")
write.csv(dfCountrybySpecies, file="CountrybySpecies.csv")


# calculate RLI for each country based on weighted species risk in that country

dfCountryRisk = data.frame(matrix(nrow=num_countries, ncol=length(all_years)+1))
colnames(dfCountryRisk) = c("Country",all_years)

for (i in 1:num_countries) {
  
  CountryRisk = filter(dfCountrybySpecies, country==df_countries$country[i])
  dfCountryRisk$Country[i] = df_countries$country[i]
  
  if (length(CountryRisk$year)>0) {
    
    years = unique(CountryRisk$year)
    
    for (j in 1:length(years)) {
      # loop through all years of record and calculate RLI
      CountryRiskbyYear = filter(CountryRisk, year==years[j])
      dfCountryRisk[[which(all_years==years[j])+1]][i] = iucn_rli(CountryRiskbyYear$risk, CountryRiskbyYear$weight)
      
    }
    
  }
  
}

save(dfCountryRisk, file="CountryRiskbyYear.Rda")
write.csv(dfCountryRisk, file="CountryRiskbyYear.csv")

dfCountryPlot = melt(dfCountryRisk, id.vars="Country", measure.vars=as.character(all_years), variable.name="year", value.name="RLI")
dfCountryPlot = na.omit(dfCountryPlot)
ggplot(dfCountryPlot, aes(x=year, y=RLI, group=Country))+geom_line()

# unique species ID list
species_id = array(0, length(all_species))

for (i in 1:length(all_species)) {
  
  species_id[i] = df_species$taxonid[which(df_species$scientific_name==all_species[i])]
  
}


iucn_rli = function(risk, weight) {
  
  # calculate Red List Index as defined by IUCN
  RLI = 1 - (sum(risk*weight) / (sum(weight)*5))
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
