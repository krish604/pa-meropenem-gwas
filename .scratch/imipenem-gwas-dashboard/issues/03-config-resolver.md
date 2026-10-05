# 03: Config resolver: science plus machine overlays

**What to build:** One place for machine-independent science and one small overlay per machine, so the laptop and the big Linux machine can be pointed at the same code without editing source, and so the two can never silently disagree about the science.

**Blocked by:** None (can start immediately)

**Status:** ready-for-agent | wired

- [ ] Scientific constants (organism, genome size expectations, breakpoints, gene and regulator tables, tool and database versions) live in exactly one machine-independent file.
- [ ] A machine overlay supplies paths, threads, memory, sample cap and per-tool availability, and nothing scientific.
- [ ] A test proves the two overlays cannot disagree about a scientific constant.
- [ ] No path, thread count, memory limit or os.cpu_count() appears anywhere in the source.

**Build state: `wired`.** Phase 2. The split is complete and every entry point goes through it: `load_config` composes science with a machine overlay, a default machine, and a refusal to silently override a self-contained file. The reference genome identity was moved OUT of the overlay after review, so two machines can no longer disagree about which genome the science uses. `os.cpu_count()` survives in three pre-existing places and is ticket 11 work.
