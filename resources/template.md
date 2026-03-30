# BCG Susceptibility Panel for NMIBC
*This test combines the results for relevant single nucleotide variants, methylation variants, and immune infiltrate prediction to reach a final score for prediction of BCG treatment success for NMIBC.*  
**Panel name** {{repo.panel_name}}<br>
**URL:** {{repo.url}}
___
### Clinical Details
**Sample Name:** {{repo.sample_name}}<br>
Bladder lesion, transurethral resection:<br>
Papillary urothelial carcinoma, high grade (grade 3/3 - WHO 1973)<br>
Noninvasive, pTa<br>
**Specimen:** TURBT biopsy  
___
## Overall Result and Interpretation
<span style="font-size:16px;">Tumour is **{{repo.word_result}}** to respond to treatment.</span><br>

*Guidance: Result is out of Highly Unlikely, Somewhat Unlikely, Somewhat Likely, Highly Likely*

**Comments:** The result is an aggregate score obtained by combining the scores of biomarkers that reach read depth threshold and fall within respective scoring parameters.

___
### Raw results
*Guidance: Total (normalised) Score is a scale of 0% – 100%. 0% will respond and 100% will not respond to BCG treatment. Most results will fall between 25% and 75%. We recommend treatment with BCG when the result is \<50% for patients who are otherwise healthy, and \<40% for all other patients.*

| **ID** | Biomarker name | Biomarker type | Result Options | Result | Score | 
|---:|---|---|---|---|---|
{% for entry in repo.all_results %}
|**{{entry[0]}}**|{{entry[1]}}|{{entry[2]}}|{{entry[3]}}|{{entry[4]}}|{{entry[5]}}|
{% endfor %}

___
### Methodology
Nanopore sequencing was used to detect single nucleotide variants, methylation status and immune status (by methylation) of the tumour sample. 

