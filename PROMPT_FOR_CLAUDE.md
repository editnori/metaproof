# Prompt for Claude Code

Paste everything below the line into Claude Code.

---

I'm doing a systematic review and meta-analysis and want to use metaproof (https://github.com/editnori/metaproof) to run it. Work through these steps in order and stop and ask me if anything fails.

1. Setup. Clone the repo to ~/Projects/metaproof and follow the "Setup on a Mac" section of its README. If something is missing (Xcode command line tools, uv, Homebrew), tell me what to run. Install R and the R engine too, since I want both engines.

2. Check the install. From the repo, run the tests (`.venv/bin/python -m pytest -q`, after installing requirements-dev.txt), the cross-check (`.venv/bin/python validation/crosscheck.py`) and the R self-check (`Rscript reference_r/tests/check_reference.R`). Tell me the counts. Anything that fails, stop and show me.

3. Read README.md, templates/HOW_TO_USE_TEMPLATES.txt and templates/protocol_template.txt so you know what each default is and why.

4. Ask me for my extraction sheet and my protocol choices: outcomes, effect measure for each, minimum number of studies, subgroups, and anything I set differently from the defaults. Don't fill in a choice I haven't made. Ask.

5. Run templates/check_extraction.py on my sheet and walk me through every error and warning. Don't fix numbers in my sheet yourself. Tell me what looks wrong and let me fix it.

6. Write the config from my choices, then run both engines on my data. Compare them. If they disagree anywhere, or if a result breaks one of my protocol rules, flag it before anything else.

7. Fill in the PRISMA 2020 flow diagram from the counts I give you. Draft the GRADE summary of findings, but the GRADE calls are mine, so show me your reasoning for each and let me decide.

8. Explain the pooled results in plain words: the effect and its unit, how sure we can be, how much the studies disagreed, and what a new study would likely show. Say what the results don't show too.

Rules:
- Use only numbers from my sheet and from metaproof's output. Never invent or estimate a number. If something's missing, tell me.
- metaproof runs on my computer, but anything you read is processed by Anthropic. Only open files my data agreements allow me to share with an AI tool.
- Don't push or publish anything without asking me.
