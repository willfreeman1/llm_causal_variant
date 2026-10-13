# Method graph

Two modes share this graph. `combined` runs the statistical model and the language model. `llm_only` skips the statistical model. Papers are retrieved only when that mode asks for them, and in the combined mode only when the two methods disagree.

```mermaid
---
config:
  flowchart:
    curve: linear
---
graph TD;
	__start__([<p>__start__</p>]):::first
	guardrail_inputs(guardrail_inputs)
	caldera(caldera)
	llm_no_papers(llm_no_papers)
	guardrail_output(guardrail_output)
	retrieve(retrieve)
	llm_with_papers(llm_with_papers)
	guardrail_output_papers(guardrail_output_papers)
	decide(decide)
	needs_review(needs_review)
	trace(trace)
	__end__([<p>__end__</p>]):::last
	__start__ --> guardrail_inputs;
	caldera --> llm_no_papers;
	decide --> needs_review;
	guardrail_inputs -.-> caldera;
	guardrail_inputs -.-> llm_no_papers;
	guardrail_inputs -.-> trace;
	guardrail_output -.-> decide;
	guardrail_output -.-> retrieve;
	guardrail_output_papers --> decide;
	llm_no_papers --> guardrail_output;
	llm_with_papers --> guardrail_output_papers;
	needs_review --> trace;
	retrieve --> llm_with_papers;
	trace --> __end__;
	classDef default fill:#f2f0ff,line-height:1.2
	classDef first fill-opacity:0
	classDef last fill:#bfb6fc

```
