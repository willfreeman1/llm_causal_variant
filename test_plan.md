# Test plan

The rules below were fixed on 9 October 2026, before any call on this test. The test was run once after that date. This page was edited afterwards for readability only. No rule changed.

A question is one trait, one lead variant, and the protein-coding genes nearby. The choosing rule is the fixed rule used on the older studies: if the language model and the statistical model name the same gene, use that gene; if they disagree, use the language model only when its sureness is at least 0.85; otherwise use the statistical model. The closest gene is the protein-coding gene nearest the lead variant.

## Which questions count

A trait is one protein measured in blood when its name matches "Circulating [name] levels", or its mapped trait matches "level of / amount of / concentration of ... in blood". Those traits are one group. The word after "Circulating" is the protein's gene symbol (`NTPROBNP` is read as `NPPB`). An Olink id in the name does not by itself name a gene.

A question is a giveaway when that symbol, or any listed gene symbol of at least 3 characters, appears as a whole word in the trait name or the mapped trait.

Repeats are the same labeled gene, leads within 500 kb, and the same trait group. One question is kept from each repeat.

- Main set: `questions/new_main.tsv`, 429 questions after repeats, with giveaways removed. 29 are lipids, 179 are blood proteins, and 221 are other traits. 150 had not already been reported before 1 September 2025. For proteins, "already reported" means the labeled gene was mapped or reported for any blood-protein trait before that date. Lipids use the chemical family. Everything else uses the same trait.
- Giveaway set: `questions/new_giveaway.tsv`, 246 questions.
- Full set: `questions/new_full.tsv`, 1,085 questions. This set is a secondary comparison. Thirteen questions that had been counted as lipids, because the mapped trait contains the word "lipoprotein", are circulating protein levels and are in the protein group. Lipids in the main set are 29.

The largest labeled gene in the main set is APOE, on 8 of 429 questions (1.9%).

The combined-system set is 38 questions. None is a giveaway. None merged with another. All 38 are inside the main set. The statistical model's picks were not compared with the labels before the run.

One neighborhood, id A00501, stays in that set of 38. Its credible set gives TREM2 a protein-changing probability of 0.41. The drop rule removes a set only when, for any gene, the summed probability of protein-changing variants is above 0.5. The check that runs on a new question uses that same rule and does not need the answer.

Two main-set representatives were members of an older question rather than the stored full-set row (M00131 and M00235). They are called as well, so the main set has an answer for the question it scores. That is 1,087 calls. The wording is the same.

## Prompts

Model `gpt-5.4`, snapshot `gpt-5.4-2026-03-05`. Temperature 0. Reasoning off. No sureness gate on the language model by itself. The system prompt is the same with and without papers. The leading newline is part of it:

```
You are an expert in biology and genetics.
Your task is to identify likely causal genes within a locus for a given GWAS phenotype based on literature evidence.
From the list, provide the likely causal gene (matching one of the given genes), confidence (0: very unsure to 1: very confident), and a brief reason (50 words or less) for your choice.
Return your response in JSON format, excluding the GWAS phenotype name and gene list in the locus. JSON keys should be 'causal_gene','confidence','reason'. 
Your response must start with '{' and end with '}'.
```

No papers:

```
Identify the causal gene.
GWAS phenotype: {trait}
Genes in locus: {gene_str}
```

With papers, on every full-set question and the two extra representatives. The 25 abstracts are the ones closest in meaning. There is no similarity cutoff. Each abstract is capped at 1,800 characters.

```
Identify the causal gene.
GWAS phenotype: {trait}
Measurement: {measurement}
Genes in locus: {gene_str}

Abstracts:
1. ({year}) {title}. {abstract}
```

Measurement sentence. A circulating protein level uses the non-lipid sentence. Lipids: the catalog trait name with "(Model 1)" or "(Model 2)" removed, then "Blood lipid measurement." Every other trait: the catalog trait name, then "Catalog term: " and the mapped trait. If those two are the same words, use the trait name only.

Paper cutoff: the earliest of the journal date and the submission date. If only one is a real date, use that one. Drop an abstract that contains the lead variant's id or its base-pair position. A paper that names the labeled gene and the trait, but not that variant, still gets through. A year with no month is kept only when that year is before the cutoff year.

On the 38-question set, the choosing rule uses 0.85. Where the two methods agree, the paper answer is ignored and the agreed gene is used. Where they disagree, the rule uses the paper answer. A missing statistical-model pick does not agree with anything. A missing answer counts as wrong.

## What is compared

- Primary, language model: with papers versus without papers, on the main set of 429.
- Primary, choosing rule: with papers versus without papers, on the 38. At 38 questions this check cannot confirm or reject a difference by itself.
- Secondary: the full 1,085, the giveaway set, the questions not previously reported, and lipids, proteins, and other traits separately. The 19 hypothyroidism questions are reported on their own and are also inside the main set.
- Every table also shows the closest gene. The statistical model alone is reported on the 38 only.

## How results are reported

- Every score is a count and a percentage, with a 95% interval from a bootstrap that redraws whole questions (2,000 draws, seed `20261009`). A second interval redraws whole labeled genes.
- Every comparison is paired: how many the new setup fixed, how many it broke, the net change with its 95% interval, and an exact two-sided McNemar p-value.
- The write-up leads with the size of the change and its range.

## Pooled comparisons

Each set alone is small. These two pools use only sets that were not used to choose the setup. Chromosomes 4, 5, 7, and 22 were used while choosing the setup and are never pooled.

1. Language model, papers versus no papers. Add the fixed and broken counts from the 18 held-out older-study chromosomes (64 fixed, 26 broken) and from the main set of 429. Report each set, the pooled counts, and an exact McNemar p-value on the pooled counts. The two sets differ: the older run added papers only where the methods disagreed, and this test adds papers on every question.
2. Choosing rule, with papers versus without. Add the 18 held-out older-study chromosomes (62 fixed, 42 broken) and the 38. Same reporting.

## Spending cap

The cap for this test is $40 in total: embeddings plus all model calls on the full set, with and without papers, including the two extra representatives. The estimate before the run was about $3.28 to embed 462,725 abstracts that were not already embedded, plus about $28.21 for the model calls, about $31.49 in all. Stop at $40 and report how many questions are left.
