function put_together()


%just to bring the values of the days when we did the experiments
round1.NVI12.days.holobmi = [190930, 191001, 191005, 191007, 191009, 191011, 191013, 191025]; 
round1.NVI13.days.holobmi = [190930, 191001, 191003, 191005, 191007, 191009, 191011, 191013, 191025]; 
round1.NVI16.days.holobmi = [190930, 191001, 191005, 191007, 191009, 191011, 191025];

round1.NVI12.days.E3 = [191004, 191026, 191028];
round1.NVI13.days.E3 = [191004, 191026, 191028];
round1.NVI16.days.E3 = [191004, 191026, 191028];

round1.NVI12.days.rr = [191006, 191010, 191014];
round1.NVI13.days.rr = [191006, 191010, 191014];
round1.NVI16.days.rr = [191006, 191010, 191014, 191024];

round1.NVI12.days.E2holo = [191008, 191012, 191027];
round1.NVI13.days.E2holo = [191008, 191012, 191027];
round1.NVI16.days.E2holo = [191008, 191012, 191027];

round1.NVI12.days.frr = [191030, 191101, 191103];
round1.NVI13.days.frr = [191030, 191101, 191103];
round1.NVI16.days.frr = [191030, 191101, 191103];

round1.NVI12.days.fhbmi = [191031, 191102, 191104];
round1.NVI13.days.fhbmi = [191031, 191102, 191104];
round1.NVI16.days.fhbmi = [191031, 191102, 191104];


round2.NVI17.days.holobmi = [191106, 191108, 191110, 191112, 191114, 191116, 191120, 191122];
round2.NVI20.days.holobmi = [191106, 191108, 191110, 191112, 191114, 191116, 191118, 191120];
round2.NVI22.days.holobmi = [191106, 191108, 191110, 191114, 191116, 191118, 191120, 191122];

round2.NVI17.days.E3 = [191109, 191115, 191121];
round2.NVI20.days.E3 = [191109, 191115, 191121];
round2.NVI22.days.E3 = [191109, 191115, 191121];

round2.NVI17.days.rr = [191111, 191117, 121123];
round2.NVI20.days.rr = [191105, 191117, 121122];
round2.NVI22.days.rr = [191105, 191111, 191117];

round2.NVI17.days.E2holo = [191107, 191113, 191119];
round2.NVI20.days.E2holo = [191107, 191113, 191119];
round2.NVI22.days.E2holo = [191107, 191113, 191119];

round2.NVI17.days.frr = [191124, 191126, 191128];
round2.NVI20.days.frr = [191124, 191126, 191128];
round2.NVI22.days.frr = [191124, 191126, 191128];

round2.NVI17.days.fhbmi = [191125, 191127, 191129];
round2.NVI20.days.fhbmi = [191125, 191127, 191130];
round2.NVI22.days.fhbmi = [191125, 191127, 191129];

%preparing the structure
animalNamesRound1 = fieldnames(round1);
animalNamesRound2 = fieldnames(round2);
animalNames = {animalNamesRound1{:}, animalNamesRound2{:}}


