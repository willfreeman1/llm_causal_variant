# Results

Each interval redraws whole questions 2,000 times (seed 20261009) and reports the middle 95 percent of those scores. Cost is the model spend for that row, divided by the number of questions.

## Practice, United Kingdom Biobank, 1,348 questions

Practice; possible training exposure. Studies in this table were public before the language model's knowledge cutoff.

| Method | Correct | % | 95% interval | Cost per question |
|---|---:|---:|---:|---:|
| Nearest gene | 636 / 1348 | 47.2 | 44.4 to 49.8 | $0 |
| Statistical model | 757 / 1348 | 56.2 | 53.4 to 58.7 | $0 |
| Language model, no papers | 798 / 1348 | 59.2 | 56.5 to 61.9 | $0.0017 |
| Combining rule, per-chromosome line | 869 / 1348 | 64.5 | 61.9 to 66.9 | $0.0017 |
| Combining rule, fixed 0.85 line | 870 / 1348 | 64.5 | 61.9 to 67.0 | $0.0017 |
| Rule with the nearest gene as the partner | 838 / 1348 | 62.2 | 59.4 to 64.7 | $0.0017 |

The per-chromosome line is the one fit on the other chromosomes. It is right on 869. A single 0.85 line on all 1,348 is right on 870. The nearest-gene partner, with its own per-chromosome line, is right on 838.

## 18 chromosomes held out of the paper-setup choice

| Method | Correct | % | 95% interval | Cost per question |
|---|---:|---:|---:|---:|
| Combining rule, no papers | 738 / 1153 | 64.0 | 61.2 to 66.8 | $0.0018 |
| Combining rule, with papers | 758 / 1153 | 65.7 | 62.9 to 68.4 | $0.014 |

Papers were requested only where the two methods disagreed. Chromosomes 4, 5, 7, and 22 are not in this table.

## New studies, main list, 429 questions

This is the primary new-study result. Studies were first published after 31 August 2025.

| Method | Correct | % | 95% interval | Cost per question |
|---|---:|---:|---:|---:|
| Nearest gene | 186 / 429 | 43.4 | 38.9 to 48.2 | $0 |
| Language model, no papers | 200 / 429 | 46.6 | 42.0 to 51.0 | $0.0018 |
| Language model, with papers | 203 / 429 | 47.3 | 42.9 to 52.0 | $0.024 |

By kind of trait, language model with no papers against the nearest gene. Counted after the run, not locked beforehand.

| Slice | Questions | Only the language model right | Only the nearest gene right | Paired p |
|---|---:|---:|---:|---:|
| All main | 429 | 95 | 81 | 0.3271 |
| Diseases and other traits | 221 | 57 | 25 | 0.0005 |
| Blood protein levels | 179 | 31 | 55 | 0.0127 |
| Not previously reported | 150 | 34 | 21 | 0.1048 |

## New studies where both methods have an input, 38 questions

Check that the setup runs; too small to conclude.

| Method | Correct | % | 95% interval | Cost per question |
|---|---:|---:|---:|---:|
| Nearest gene | 14 / 38 | 36.8 | 21.1 to 52.6 | $0 |
| Language model, no papers | 23 / 38 | 60.5 | 44.7 to 76.3 | $0.0019 |
| Language model, with papers | 24 / 38 | 63.2 | 47.4 to 79.0 | $0.024 |
| Statistical model alone | 12 / 38 | 31.6 | 18.4 to 47.4 | $0 |
| Combining rule, no papers | 24 / 38 | 63.2 | 47.4 to 79.0 | $0.0019 |
| Combining rule, with papers | 24 / 38 | 63.2 | 47.4 to 79.0 | $0.026 |

## Earlier checks

List order and fame, on 80 neighborhoods asked four times each. The top gene was the same on every order for 78 of 80. Putting the true gene in the first third of the list was right on 93.8 percent of those trials, and in the last third on 92.5 percent. When the answer was wrong, the chosen gene was more studied than the truth 50 percent of the time, against 53 percent for a random other gene in the same list. Spend $0.601.

Memorization probe. Exact copy of a real Open Targets row: 0 of 25. The same question on the iris flower file, which models often have memorized: 2 of 15. Spend $0.106.

## Guardrails

Answers outside the candidate list were stored as blank and counted wrong: 77 with no papers, 5 with papers. All of them were blood-protein questions, and in every one the model named the measured protein.
Abstracts containing the lead variant's id or position were dropped when the pools were built: 506 on the new-study questions and 0 on the held-out chromosomes. Retrieval is LanceDB vector search with gene and date prefilters. On neighborhood A00407 the same 25 papers come back as in the locked run, with two swapped, because those scores are equal to 7 decimals. The locked run did not save a separate count of how many abstracts the date alone removed.
Statistical-model sets dropped for a protein-changing leak before the run: 7 (coding probability over 0.5 on 2, the protein-changing lead inside the set on 6).
Questions flagged for a person to review, without changing the pick: 106 of 1,348 on the practice set (27 correct, 25.47 percent) against 67.79 percent on the rest. On the new main list with no papers, 118 flagged (18.64 percent correct) against 57.23 percent on the rest.
One kept input, TREM2, has a protein-changing probability of 0.41. The locked run kept it because the drop line was 0.5, and both methods named TREM2, so it is one of the 24. A new run uses that same line and does not look at the answer, so it keeps the question too.

Total API spend recorded for the project: $72.634.

On the new main list, adding papers changed the language model by 3 questions, and redrawing the questions puts that change between a loss of 11 and a gain of 18. The same model, with no papers, beats the nearest gene on diseases and other traits and loses to it on blood protein levels. The combining rule's paper change on the 38 new questions that have both inputs was zero. On the 18 held-out chromosomes, papers changed the rule by 20 questions, and those studies may have been in the training data. Counting the same protein gene once, instead of once per protein trait, is what keeps the new result from being dominated by one repeated fact.

## Diagnostic run after the locked test, on questions already seen

One sentence was added to the instructions, telling the model that the measured protein's gene may be elsewhere and that the question is which listed gene changes the measurement. No papers. The locked scores above are unchanged. Cost $2.148.
Blood protein levels, main set: 73 to 83 of 179. Fixed 14, broken 4, net +10. Redrawing the questions puts that change between +2 and +18.
Diseases and other traits: 110 to 112 of 221, net +2. Lipids: 17 of 29 either way. Main set: 200 to 212 of 429, net +12. Giveaway questions, where the protein's own gene is the right answer: 228 of 246 either way.
Answers outside the list fell from 77 to 26 overall, and from 23 to 12 on main-set protein questions. On those protein questions the model was the only correct method on 40 and the nearest gene was the only correct method on 54. The nearest gene, at 97 of 179, is still ahead of the model, at 83.
On the main set, the median sureness moved from 0.93 to 0.86. Answers at or above 0.85, the line the combining rule uses, moved from 251 to 224 of 429. The combining rule was not re-run.
