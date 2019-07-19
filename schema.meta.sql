CREATE TABLE subject
  ( subject_id PRIMARY KEY
  , ibd_diagnosis_at_age
  , age
  , site
  , history_other_abdominal_surgery
  , education
  , occupation
  , history_tonsillectomy
  , ibd_diagnosis
  , history_childhood_farm
  , history_childhood_daycare
  , history_childhood_environmental_tobacco_smoke
  , history_birth_premature
  , history_birth_hospital
  , history_birth_csection
  , history_childhood_breastfed
  , history_childhood_antibiotics_pre_age_one
  , history_childhood_antibiotics_pre_age_five
  , history_childhood_pets
  , race
  , sex
  , status_smoker
  , status_smoker_years
  , history_smoker_age_at_start
  , status_smoker_number_per_day
  , baseline_height
  , baseline_weight
  , history_smoker
  );

CREATE TABLE visit
  ( visit_id PRIMARY KEY
  , subject_id REFERENCES subject(subject_id)
  , week_number
  , visit_date
  , visit_number
  , diet_sugary_beverage
  , diet_artificial_sweetener_beverage
  , diet_fruit_juice
  , diet_water
  , diet_alcohol
  , diet_live_bacteria
  , diet_dairy
  , diet_probiotic
  , diet_whole_fruits
  , diet_whole_vegetables
  , diet_beans
  , diet_whole_grains
  , diet_starch
  , diet_eggs
  , diet_meats_processed
  , diet_meats_red
  , diet_meats_white
  , diet_meats_shellfish
  , diet_meats_fish
  , diet_sweets
  , status_antibiotics
  , status_chemotherapy
  , status_immunosuppressants
  , status_recent_colonoscopy
  , status_oral_contrast
  , status_diarrhea
  , status_recent_hospitalization
  , status_ever_had_bowel_surgery
  , diet_tea_or_coffee
  , status_general_wellbeing
  , status_sccai
  , status_defecation_urgency
  , status_blood_in_stool
  , status_abdominal_pain
  , status_number_soft_stools
  , status_abdominal_mass
  , status_arthralgia
  , status_hbi
  , status_crp
  , status_esr
  , status_weight
  );

CREATE TABLE stool
  ( stool_id PRIMARY KEY
  , visit_id REFERENCES visit(visit_id)
  , fecal_calprotectin
  );

CREATE TABLE preparation
  ( preparation_id PRIMARY KEY
  , project_name
  , library_type
  , stool_id REFERENCES stool(stool_id)
  );

CREATE TABLE library
  ( library_id PRIMARY KEY
  , external_id
  , library_type
  , number_of_lanes_aggregated
  , sequenced_reads
  , preparation_id REFERENCES preparation(preparation_id)
  );

CREATE TABLE library_group
  ( library_group
  , library_id
  );

CREATE VIEW mgen_library AS
SELECT * FROM library WHERE library_type = 'metagenomics'
;
