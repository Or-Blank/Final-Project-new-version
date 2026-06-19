from __future__ import annotations

import os
import sys
import threading
import warnings
from typing import Optional

try:
    import tkinter as tk
    from tkinter import filedialog, messagebox, scrolledtext, ttk
    TK_AVAILABLE = True
except ImportError:
    TK_AVAILABLE = False

from .loader import BCRTreeLoader
from .pipeline import run as run_pipeline


def _info_btn(parent: tk.Widget, popup_title: str, popup_text: str) -> tk.Label:
    """Return a clickable  ℹ️  label that opens a messagebox popup."""
    lbl = tk.Label(parent, text=" ℹ️ ", font=("Helvetica", 10),
                   cursor="hand2", fg="#2c7bb6")
    lbl.bind("<Button-1>", lambda _e: messagebox.showinfo(popup_title, popup_text))
    return lbl


def _open_clone_picker(parent: tk.Tk, clone_sizes: dict,
                       current: str, on_select) -> None:
    """Modal scrollable list of clones with >= 2 cells, sorted by size."""

    # Only show clones with at least 2 observed cells (singletons excluded)
    eligible = {cid: n for cid, n in clone_sizes.items() if n >= 2}

    win = tk.Toplevel(parent)
    win.title("Select Clone")
    win.geometry("440x540")
    win.resizable(True, True)
    win.grab_set()  # make modal

    tk.Label(
        win,
        text=(f"Showing {len(eligible)} clones with ≥ 2 cells  "
              f"(out of {len(clone_sizes)} total in file)."),
        font=("Helvetica", 9), fg="#555",
    ).pack(pady=(10, 4), padx=10, anchor="w")

    # ── search / filter row ───────────────────────────────────────────────
    search_var = tk.StringVar()
    sf = tk.Frame(win)
    sf.pack(fill="x", padx=10, pady=(0, 4))
    tk.Label(sf, text="🔍  Filter:").pack(side="left")
    search_entry = ttk.Entry(sf, textvariable=search_var, width=30)
    search_entry.pack(side="left", padx=6)
    search_entry.focus_set()

    # ── listbox + scrollbar ───────────────────────────────────────────────
    lf = tk.Frame(win)
    lf.pack(fill="both", expand=True, padx=10)
    sb = ttk.Scrollbar(lf, orient="vertical")
    lb = tk.Listbox(lf, yscrollcommand=sb.set, font=("Courier", 9),
                    selectmode="single", activestyle="dotbox", height=22)
    sb.config(command=lb.yview)
    sb.pack(side="right", fill="y")
    lb.pack(side="left", fill="both", expand=True)

    sorted_clones = sorted(eligible.items(), key=lambda x: x[1], reverse=True)

    def _populate(flt: str = "") -> None:
        lb.delete(0, "end")
        for cid, n in sorted_clones:
            if flt.lower() in cid.lower():
                lb.insert("end", f"{cid:<38}  {n:>5} cells")
        # Pre-select the previously chosen clone if it is still visible
        for i in range(lb.size()):
            if lb.get(i).startswith(current):
                lb.selection_set(i)
                lb.see(i)
                break

    _populate()
    search_var.trace_add("write", lambda *_: _populate(search_var.get()))

    info_var = tk.StringVar(value=f"{len(eligible)} eligible clones")
    tk.Label(win, textvariable=info_var,
             font=("Helvetica", 8), fg="#888").pack(pady=(2, 0))

    def _on_sel(*_args) -> None:
        sel = lb.curselection()
        if sel:
            cid = lb.get(sel[0]).split()[0]
            info_var.set(f"Selected: {cid}  ({eligible.get(cid, '?')} cells)")

    lb.bind("<<ListboxSelect>>", _on_sel)

    # ── buttons ───────────────────────────────────────────────────────────
    bf = tk.Frame(win)
    bf.pack(fill="x", padx=10, pady=8)

    def _confirm() -> None:
        sel = lb.curselection()
        if sel:
            on_select(lb.get(sel[0]).split()[0])
        win.destroy()

    def _all_clones() -> None:
        on_select(None)
        win.destroy()

    ttk.Button(bf, text="✓  Use Selected Clone",
               command=_confirm).pack(side="left", ipadx=8)
    ttk.Button(bf, text="⊞  All Clones",
               command=_all_clones).pack(side="left", padx=8)
    ttk.Button(bf, text="Cancel",
               command=win.destroy).pack(side="right")
    lb.bind("<Double-Button-1>", lambda _e: _confirm())


# ── Popup text constants ──────────────────────────────────────────────────────

_REQUIRED_COLS_TEXT = (
    "ALWAYS REQUIRED\n"
    "─────────────────────────────────────────\n"
    "  clone_id\n"
    "      Clonal family identifier.\n"
    "      Aliases: clonotype, clonotype_id, lineage_id\n\n"
    "  VDJ_sequence_H\n"
    "      Heavy-chain VDJ nucleotide sequence.\n"
    "      Aliases: sequence_alignment_H, sequence_alignment,\n"
    "               heavy_sequence\n\n"
    "STRONGLY RECOMMENDED (warning if missing)\n"
    "─────────────────────────────────────────\n"
    "  germline_alignment_d_mask\n"
    "      IMGT-masked germline reference sequence.\n"
    "      Without it a consensus is estimated from observed\n"
    "      sequences (less accurate).\n"
    "      Aliases: germline_alignment, uca_sequence_H\n\n"
    "OPTIONAL  (tree still builds; metadata shows 'Unknown')\n"
    "─────────────────────────────────────────\n"
    "  VDJ_sequence_L\n"
    "      Light-chain sequence. H+L concatenated for better\n"
    "      accuracy when present.\n"
    "      Alias: sequence_alignment_L, light_sequence\n\n"
    "  c_call\n"
    "      Isotype call (e.g. IGHG1, IGHM).\n"
    "      Used for node colours and CSR refinement.\n"
    "      Aliases: isotype, ig_class\n\n"
    "  sample_id\n"
    "      Tissue or sample of origin.\n"
    "      Aliases: Sample, tissue, organ, tissue_origin\n\n"
    "  cluster_annotated\n"
    "      Cell-type / GC-subset label.\n"
    "      Aliases: annot_clusters, cell_type, leiden\n\n"
    "  cell_id\n"
    "      Unique cell identifier / barcode.\n"
    "      Falls back to row number if absent.\n"
    "      Aliases: barcode, sequence_id\n\n"
    "ORDER OF COLUMNS\n"
    "─────────────────────────────────────────\n"
    "Column order DOES NOT matter.\n"
    "IMGT-gapped sequences (dots/dashes) are handled automatically."
)

_MAX_CLONES_TEXT = (
    "WHAT IT DOES?\n\n"
    "Limits how many clonal families are processed in one run.\n\n"
    "IS IT NECESSARY?\n\n"
    "No. Leave blank to process every eligible clone in the file.\n\n"
    "WHEN TO USE IT:\n\n"
    "  • Quick test: set to 5 to check the pipeline works on\n"
    "    your file before committing to a full run.\n"
    "  • Large files (1000+ clones) may take several minutes;\n"
    "    a small limit lets you preview results immediately.\n\n"
    "EXAMPLE:\n\n"
    "  5     →  only the first 5 eligible clones are processed.\n"
    "  blank →  all eligible clones."
)

_COLLAPSE_TEXT = (
    "WHAT IT DOES?\n\n"
    "B cell sequences are short (~300–700 bp) and differ from each\n"
    "other by only a few mutations. This means the tree algorithm\n"
    "often cannot decide whether two cells split before or after a\n"
    "mutation — it produces a near-zero branch length between them.\n\n"
    "This setting collapses those ambiguous branches into a single\n"
    "fork (polytomy), which is biologically more honest than forcing a fake ordering.\n\n"
    "IS IT NECESSARY?\n\n"
    "No. The default (1e-6) is appropriate for virtually all BCR\n"
    "datasets and should not be changed unless you have a specific reason.\n\n"
    "WHEN TO CHANGE IT:\n\n"
    "  • Increase (e.g. 0.001) if trees look too bushy and you want\n"
    "    to merge more branches.\n"
    "  • Decrease (e.g. 0) to keep every branch even if the\n"
    "    resolution is statistically uncertain.\n\n"
    "EXAMPLE VALUES:\n\n"
    "  1e-6  →  default, removes floating-point noise only\n"
    "  0.001 →  aggressive collapsing, simpler tree\n"
    "  0     →  no collapsing at all"
)

_NNI_TEXT = (
    "WHAT IT DOES?\n\n"
    "After building the initial tree from sequence distances, this\n"
    "step checks whether any parent→child edges imply a biologically impossible isotype switch.\n\n"
    "B cells can only switch isotype in ONE direction:\n\n"
    "IgM → IgD → IgG3 → IgG1 → IgA1 → IgG2 → IgG4 → IgE → IgA2\n\n"
    "This is called Class Switch Recombination (CSR) and is\n"
    "irreversible. If the sequence-based tree places an IgG cell as\n"
    "the parent of an IgM cell, that is a mistake.\n\n"
    "The refinement performs local swaps (NNI = Nearest-Neighbour Interchange) to fix those violations, using the sequence signal as a tiebreaker so the tree is not changed more than necessary.\n\n"
    "IS IT NECESSARY?\n\n"
    "Recommended ON. The improvement is small when the sequence signal is strong, but useful when many cells share the same isotype or when sequences are very similar.\n\n"
    "WHEN TO TURN IT OFF:\n\n"
    "  • Your file has no c_call / isotype column.\n"
    "  • You want the pure sequence-based tree for comparison.\n"
    "  • Speed: for very large clones (500+ cells) this step\n"
    "    adds a few seconds."
)


# ── Main GUI ──────────────────────────────────────────────────────────────────

def launch_gui() -> None:
    """Open the BCR Lineage Tracer GUI.  Blocks until the window is closed."""

    if not TK_AVAILABLE:
        print(
            "tkinter is not available.\n"
            "  Windows / macOS : reinstall Python from python.org\n"
            "  Ubuntu / Debian : sudo apt install python3-tk\n"
            "  Fedora          : sudo dnf install python3-tkinter"
        )
        sys.exit(1)

    # ── mutable state ─────────────────────────────────────────────────────
    _clone_sizes: dict        = {}
    _selected_clone: Optional[str] = None

    # ── root window ───────────────────────────────────────────────────────
    root = tk.Tk()
    root.title("BCR Lineage Tracer")
    root.geometry("840x800")
    root.resizable(True, True)

    style = ttk.Style()
    style.theme_use("clam")
    style.configure("Section.TLabelframe.Label",
                    font=("Helvetica", 10, "bold"))
    style.configure("Run.TButton", font=("Helvetica", 11, "bold"))
    style.configure("Toggle.TButton", font=("Helvetica", 9))

    # ── header bar ────────────────────────────────────────────────────────
    hdr = tk.Frame(root, bg="#1a3a5c", height=58)
    hdr.pack(fill="x")
    tk.Label(hdr, text="BCR Lineage Tracer 🧬🌳",
             font=("Helvetica", 18, "bold"),
             fg="white", bg="#1a3a5c").pack(side="left", padx=16, pady=10)
    tk.Label(hdr, text="B cell clonal lineage tree maker",
             font=("Helvetica", 10), fg="#aecbe8",
             bg="#1a3a5c").pack(side="left", pady=18)

    main = tk.Frame(root, bg="#f5f5f5")
    main.pack(fill="both", expand=True, padx=12, pady=8)

    # ═══════════════════════════════════════════════════════════════════════
    # SECTION 1 — File input
    # ═══════════════════════════════════════════════════════════════════════
    inp_frame = ttk.LabelFrame(main, text="Input / Output",
                               style="Section.TLabelframe", padding=8)
    inp_frame.pack(fill="x", pady=4)
    inp_frame.columnconfigure(1, weight=1)

    inp_var = tk.StringVar()

    tk.Label(inp_frame, text="Excel file (.xlsx):").grid(
        row=0, column=0, sticky="w")
    ttk.Entry(inp_frame, textvariable=inp_var, width=50).grid(
        row=0, column=1, padx=6, sticky="ew")
    _info_btn(inp_frame, "Required Excel Columns",
              _REQUIRED_COLS_TEXT).grid(row=0, column=2, padx=(4, 0))
    ttk.Button(inp_frame, text="Browse…",
               command=lambda: _browse_in()).grid(row=0, column=3,
                                                  padx=(4, 0))

    out_var = tk.StringVar(value="bcr_lineage_output")
    tk.Label(inp_frame, text="Output directory:").grid(
        row=1, column=0, sticky="w", pady=(6, 0))
    ttk.Entry(inp_frame, textvariable=out_var, width=50).grid(
        row=1, column=1, padx=6, sticky="ew", pady=(6, 0))
    ttk.Button(inp_frame, text="Browse…",
               command=lambda: _browse_out()).grid(row=1, column=3,
                                                   padx=(4, 0),
                                                   pady=(6, 0))

    def _browse_in() -> None:
        p = filedialog.askopenfilename(
            filetypes=[("Excel files", "*.xlsx"), ("All files", "*.*")])
        if p:
            inp_var.set(p)
            _load_clone_list(p)

    def _browse_out() -> None:
        p = filedialog.askdirectory()
        if p:
            out_var.set(p)

    # ═══════════════════════════════════════════════════════════════════════
    # SECTION 2 — Clone selection
    # ═══════════════════════════════════════════════════════════════════════
    clone_frame = ttk.LabelFrame(main, text="Clone Selection",
                                 style="Section.TLabelframe", padding=8)
    clone_frame.pack(fill="x", pady=4)
    clone_frame.columnconfigure(1, weight=1)

    clone_label_var = tk.StringVar(
        value="Load a file first, then click 'Select Clone…'")
    tk.Label(clone_frame, textvariable=clone_label_var,
             font=("Helvetica", 9), fg="#444",
             anchor="w").grid(row=0, column=0, columnspan=2,
                              sticky="ew", padx=(0, 6))

    def _load_clone_list(filepath: str) -> None:
        nonlocal _clone_sizes, _selected_clone
        try:
            with warnings.catch_warnings():
                warnings.simplefilter("ignore")
                loader = BCRTreeLoader(filepath).load()
            clones = loader.get_clones()
            _clone_sizes = {
                cid: sum(1 for r in recs if not r.is_germline)
                for cid, recs in clones.items()
            }
            _selected_clone = None
            n_total    = len(_clone_sizes)
            n_eligible = sum(1 for n in _clone_sizes.values() if n >= 2)
            clone_label_var.set(
                f"File loaded  ·  {n_total} clones found  "
                f"·  {n_eligible} eligible (≥ 2 cells)  "
                f"·  Currently: ALL eligible clones"
            )
            pick_btn.config(state="normal")
        except Exception as exc:     # noqa: BLE001
            clone_label_var.set(f"⚠  Could not read file: {exc}")
            _clone_sizes = {}
            pick_btn.config(state="disabled")

    def _pick_clone() -> None:
        if not _clone_sizes:
            messagebox.showinfo("No file loaded",
                                "Please load a file first.")
            return

        def _on_pick(cid: Optional[str]) -> None:
            nonlocal _selected_clone
            _selected_clone = cid
            if cid is None:
                n_el = sum(1 for n in _clone_sizes.values() if n >= 2)
                clone_label_var.set(
                    f"{len(_clone_sizes)} clones  ·  "
                    f"Processing: ALL {n_el} eligible clones")
            else:
                n = _clone_sizes.get(cid, "?")
                clone_label_var.set(
                    f"{len(_clone_sizes)} clones  ·  "
                    f"Selected: {cid}  ({n} cells)")

        _open_clone_picker(root, _clone_sizes,
                           current=_selected_clone or "",
                           on_select=_on_pick)

    def _clear_clone() -> None:
        nonlocal _selected_clone
        _selected_clone = None
        n     = len(_clone_sizes)
        n_el  = sum(1 for v in _clone_sizes.values() if v >= 2)
        clone_label_var.set(
            (f"{n} clones  ·  Processing: ALL {n_el} eligible clones"
             if n else
             "Load a file first, then click 'Select Clone…'"))

    pick_btn = ttk.Button(clone_frame, text="🔍  Select Clone…",
                          command=_pick_clone, state="disabled")
    pick_btn.grid(row=0, column=2, padx=(8, 0))

    ttk.Button(clone_frame, text="✕  Clear (use all)",
               command=_clear_clone).grid(row=0, column=3, padx=(4, 0))

    # ═══════════════════════════════════════════════════════════════════════
    # SECTION 3 — Advanced options  (collapsible)
    # ═══════════════════════════════════════════════════════════════════════
    adv_outer   = tk.Frame(main, bg="#f5f5f5")
    adv_outer.pack(fill="x", pady=2)
    adv_visible = tk.BooleanVar(value=False)
    adv_inner   = ttk.LabelFrame(adv_outer, text="Advanced Options",
                                  style="Section.TLabelframe", padding=10)

    def _toggle_adv() -> None:
        if adv_visible.get():
            adv_inner.pack_forget()
            adv_visible.set(False)
            adv_toggle_btn.config(
                text="▸  Advanced Options  (click to expand)")
        else:
            adv_inner.pack(fill="x")
            adv_visible.set(True)
            adv_toggle_btn.config(
                text="▾  Advanced Options  (click to collapse)")

    adv_toggle_btn = ttk.Button(
        adv_outer,
        text="▸  Advanced Options  (click to expand)",
        command=_toggle_adv,
        style="Toggle.TButton",
    )
    adv_toggle_btn.pack(anchor="w")

    # ── helper: one consistent option row ─────────────────────────────────
    # Each row has:  [label]  [entry/checkbox]  [ℹ️]  [optional tip]
    # This keeps all three options visually identical.

    def _option_row(parent: tk.Widget, label: str,
                    info_title: str, info_text: str) -> tk.Frame:
        """Return a frame containing a bold label and ℹ️ button."""
        row = tk.Frame(parent, bg="#f5f5f5")
        tk.Label(row, text=label, font=("Helvetica", 9, "bold"),
                 anchor="w", width=38, bg="#f5f5f5").pack(side="left")
        _info_btn(row, info_title, info_text).pack(side="left", padx=2)
        return row

    # ── Max clones ────────────────────────────────────────────────────────
    max_var = tk.StringVar(value="")

    max_row = _option_row(adv_inner,
                          "Max clones   (blank = process all)",
                          "Max Clones", _MAX_CLONES_TEXT)
    max_row.pack(fill="x", pady=(0, 6))
    ttk.Entry(max_row, textvariable=max_var,
              width=8).pack(side="left", padx=6)
    tk.Label(max_row, text="💡 Set to 5 for a quick test.",
             font=("Helvetica", 8), fg="#888",
             bg="#f5f5f5").pack(side="left", padx=8)

    ttk.Separator(adv_inner, orient="horizontal").pack(fill="x", pady=4)

    # ── Collapse threshold ────────────────────────────────────────────────
    thr_var = tk.StringVar(value="1e-6")

    thr_row = _option_row(adv_inner,
                          "Collapse threshold   (default 1e-6)",
                          "Collapse Threshold", _COLLAPSE_TEXT)
    thr_row.pack(fill="x", pady=(0, 6))
    ttk.Entry(thr_row, textvariable=thr_var,
              width=10).pack(side="left", padx=6)
    tk.Label(thr_row,
             text="💡 Leave at default unless trees look wrong.",
             font=("Helvetica", 8), fg="#888",
             bg="#f5f5f5").pack(side="left", padx=8)

    ttk.Separator(adv_inner, orient="horizontal").pack(fill="x", pady=4)

    # ── Isotype-aware NNI ─────────────────────────────────────────────────
    # Custom checkbox that shows ✓ / ☐ instead of the OS default symbol
    # (which appears as an X on many Windows themes — fix 5).
    refine_var  = tk.BooleanVar(value=True)
    check_text  = tk.StringVar(
        value="✓  Isotype-aware NNI refinement  (recommended ON)")

    def _on_refine_toggle(*_args) -> None:
        if refine_var.get():
            check_text.set("✓  Isotype-aware NNI refinement  (recommended ON)")
        else:
            check_text.set("☐  Isotype-aware NNI refinement  (recommended ON)")

    refine_var.trace_add("write", _on_refine_toggle)

    ref_row = tk.Frame(adv_inner, bg="#f5f5f5")
    ref_row.pack(fill="x", pady=(0, 4))

    tk.Checkbutton(
        ref_row,
        textvariable=check_text,
        variable=refine_var,
        font=("Helvetica", 9),
        bg="#f5f5f5",
        activebackground="#f5f5f5",
        indicatoron=False,   # hide the OS-drawn indicator box entirely
        relief="flat",
        bd=0,
        padx=0,
    ).pack(side="left")

    _info_btn(ref_row,
              "Isotype-aware NNI Refinement",
              _NNI_TEXT).pack(side="left", padx=4)

    tk.Label(adv_inner,
             text="Click ℹ️ next to each option for a full explanation.",
             font=("Helvetica", 8), fg="#aaa",
             bg="#f5f5f5").pack(anchor="w", pady=(6, 0))

    # ═══════════════════════════════════════════════════════════════════════
    # SECTION 4 — Log
    # ═══════════════════════════════════════════════════════════════════════
    log_frame = ttk.LabelFrame(main, text="Log",
                               style="Section.TLabelframe", padding=4)
    log_frame.pack(fill="both", expand=True, pady=4)

    log_box = scrolledtext.ScrolledText(
        log_frame, height=12, font=("Courier", 9),
        bg="#1e1e2e", fg="#cdd6f4", insertbackground="white",
        state="disabled",
    )
    log_box.pack(fill="both", expand=True)

    def log_msg(msg: str) -> None:
        log_box.config(state="normal")
        log_box.insert("end", msg + "\n")
        log_box.see("end")
        log_box.config(state="disabled")
        root.update_idletasks()

    prog = ttk.Progressbar(main, mode="indeterminate")
    prog.pack(fill="x", pady=(0, 4))

    # ═══════════════════════════════════════════════════════════════════════
    # SECTION 5 — Run buttons
    # ═══════════════════════════════════════════════════════════════════════
    btn_frame = tk.Frame(main, bg="#f5f5f5")
    btn_frame.pack(fill="x")

    run_btn  = ttk.Button(btn_frame, text="▶  Run Pipeline",
                          style="Run.TButton")
    run_btn.pack(side="left", ipadx=16, ipady=4)

    open_btn = ttk.Button(btn_frame, text="📂  Open Output Folder",
                          state="disabled")
    open_btn.pack(side="left", padx=8)

    # ── run logic (background thread keeps UI responsive) ─────────────────
    def _do_run() -> None:
        inp = inp_var.get().strip()
        if not inp:
            messagebox.showerror("No file selected",
                                 "Please select an Excel (.xlsx) file.")
            return
        if not os.path.exists(inp):
            messagebox.showerror("File not found",
                                 f"Cannot find:\n{inp}")
            return

        out = out_var.get().strip() or "bcr_lineage_output"
        cid = _selected_clone

        try:
            thr = float(thr_var.get())
        except ValueError:
            thr = 1e-6
        try:
            mc = int(max_var.get()) if max_var.get().strip() else None
        except ValueError:
            mc = None

        log_box.config(state="normal")
        log_box.delete("1.0", "end")
        log_box.config(state="disabled")

        run_btn.config(state="disabled")
        open_btn.config(state="disabled")
        prog.start(10)

        log_msg(f"Input  : {inp}")
        log_msg(f"Output : {out}")
        log_msg(f"Clone  : {cid if cid else 'ALL eligible'}")
        log_msg("─" * 60)

        def _worker() -> None:
            try:
                ok, skip, fail, _ = run_pipeline(
                    inp, out,
                    clone_id=cid,
                    collapse_threshold=thr,
                    refine_isotypes=refine_var.get(),
                    max_clones=mc,
                    progress_cb=log_msg,
                )
                log_msg("─" * 60)
                log_msg(f"✅  Complete — {ok} trees built, "
                        f"{skip} skipped (<2 cells), {fail} failed.")

                def _open_folder() -> None:
                    if sys.platform == "win32":
                        os.startfile(out)          # type: ignore[attr-defined]
                    elif sys.platform == "darwin":
                        os.system(f'open "{out}"')
                    else:
                        os.system(f'xdg-open "{out}"')

                open_btn.config(state="normal", command=_open_folder)

            except Exception as exc:               # noqa: BLE001
                log_msg(f"❌  Error: {exc}")
            finally:
                prog.stop()
                run_btn.config(state="normal")

        threading.Thread(target=_worker, daemon=True).start()

    run_btn.config(command=_do_run)
    root.mainloop()