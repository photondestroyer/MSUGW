"""Small layout engine (reportlab) for the plan PDF: headings, paragraphs, tables, step cards, figures."""
import os
import re

from reportlab.lib import colors
from reportlab.lib.enums import TA_LEFT
from reportlab.lib.pagesizes import A4
from reportlab.lib.styles import ParagraphStyle
from reportlab.lib.units import mm
from reportlab.pdfbase import pdfmetrics
from reportlab.pdfbase.ttfonts import TTFont
from reportlab.platypus import (SimpleDocTemplate, Paragraph, Spacer, Table, TableStyle, Image, PageBreak,
                                KeepTogether, CondPageBreak)

FONTS = "C:/Windows/Fonts"
pdfmetrics.registerFont(TTFont("A", f"{FONTS}/arial.ttf"))
pdfmetrics.registerFont(TTFont("A-B", f"{FONTS}/arialbd.ttf"))
pdfmetrics.registerFont(TTFont("A-I", f"{FONTS}/ariali.ttf"))
pdfmetrics.registerFont(TTFont("A-BI", f"{FONTS}/arialbi.ttf"))
pdfmetrics.registerFontFamily("A", normal="A", bold="A-B", italic="A-I", boldItalic="A-BI")

INK, INK2, MUTED = colors.HexColor("#0b0b0b"), colors.HexColor("#3f3e3b"), colors.HexColor("#6b6964")
NAVY, BLUE, PALE, WARM, RULE = (colors.HexColor("#1c4f8f"), colors.HexColor("#2a78d6"), colors.HexColor("#e8f1fc"),
                                colors.HexColor("#f4f2ec"), colors.HexColor("#c9c7bf"))
PAGE_W, PAGE_H = A4
LM = RM = 17 * mm
TM, BM = 17 * mm, 16 * mm
WIDTH = PAGE_W - LM - RM

ST = dict(
    title=ParagraphStyle("title", fontName="A-B", fontSize=19, leading=23, textColor=NAVY, spaceAfter=3),
    sub=ParagraphStyle("sub", fontName="A", fontSize=11.5, leading=15, textColor=INK2, spaceAfter=4),
    meta=ParagraphStyle("meta", fontName="A", fontSize=8.2, leading=11, textColor=MUTED, spaceAfter=7),
    h1=ParagraphStyle("h1", fontName="A-B", fontSize=13.2, leading=16.5, textColor=NAVY, spaceBefore=11, spaceAfter=4.5),
    h2=ParagraphStyle("h2", fontName="A-B", fontSize=10.6, leading=13.5, textColor=INK, spaceBefore=8, spaceAfter=3),
    body=ParagraphStyle("body", fontName="A", fontSize=9.1, leading=12.3, textColor=INK, spaceAfter=4.5, alignment=TA_LEFT),
    bullet=ParagraphStyle("bullet", fontName="A", fontSize=9.1, leading=12.2, textColor=INK, leftIndent=11, bulletIndent=2,
                          spaceAfter=2.2),
    cell=ParagraphStyle("cell", fontName="A", fontSize=7.9, leading=9.9, textColor=INK),
    cellh=ParagraphStyle("cellh", fontName="A-B", fontSize=7.9, leading=9.9, textColor=colors.white),
    card=ParagraphStyle("card", fontName="A", fontSize=8.35, leading=10.7, textColor=INK),
    cardlab=ParagraphStyle("cardlab", fontName="A-B", fontSize=7.4, leading=10.7, textColor=MUTED),
    cardhead=ParagraphStyle("cardhead", fontName="A-B", fontSize=9.6, leading=12, textColor=colors.white),
    cap=ParagraphStyle("cap", fontName="A-I", fontSize=8.0, leading=10.3, textColor=INK2, spaceBefore=2, spaceAfter=7),
    eq=ParagraphStyle("eq", fontName="A", fontSize=9.6, leading=13, textColor=INK, leftIndent=16, spaceBefore=1, spaceAfter=5),
    ref=ParagraphStyle("ref", fontName="A", fontSize=7.7, leading=9.8, textColor=INK, leftIndent=10, firstLineIndent=-10,
                       spaceAfter=2.0),
)
USED_TEXT = []


def M(text):
    """Plain text with light markup -> reportlab XML. **bold**  //italic//  ^{sup}  _{sub}"""
    USED_TEXT.append(text)
    t = text.replace("&", "&amp;").replace("<", "&lt;").replace(">", "&gt;")
    t = re.sub(r"\*\*(.+?)\*\*", r"<b>\1</b>", t)
    t = re.sub(r"//(.+?)//", r"<i>\1</i>", t)
    t = re.sub(r"\^\{(.+?)\}", r"<super>\1</super>", t)
    t = re.sub(r"_\{(.+?)\}", r"<sub>\1</sub>", t)
    return t.replace("\n", "<br/>")


def P(text, style="body"):
    return Paragraph(M(text), ST[style])


def H1(text):
    return [CondPageBreak(40 * mm), Paragraph(M(text), ST["h1"])]


def H2(text):
    return [CondPageBreak(26 * mm), Paragraph(M(text), ST["h2"])]


def BUL(items):
    return [Paragraph(M(i), ST["bullet"], bulletText="\u2022") for i in items]


def SP(h=3):
    return Spacer(1, h * mm)


def TABLE(header, rows, widths, caption=None, align_center=(), font=None, keep=False):
    """widths in mm (scaled to the text width). Cells are marked-up strings."""
    scale = WIDTH / (sum(widths) * mm)
    cw = [w * mm * scale for w in widths]
    cs = ST["cell"] if font is None else ParagraphStyle("c2", parent=ST["cell"], fontSize=font, leading=font + 2.0)
    hs = ST["cellh"] if font is None else ParagraphStyle("h2c", parent=ST["cellh"], fontSize=font, leading=font + 2.0)
    data = [[Paragraph(M(h), hs) for h in header]] + [[Paragraph(M(str(c)), cs) for c in r] for r in rows]
    t = Table(data, colWidths=cw, repeatRows=1)
    style = [("BACKGROUND", (0, 0), (-1, 0), NAVY), ("VALIGN", (0, 0), (-1, -1), "TOP"),
             ("LEFTPADDING", (0, 0), (-1, -1), 3.2), ("RIGHTPADDING", (0, 0), (-1, -1), 3.2),
             ("TOPPADDING", (0, 0), (-1, -1), 2.2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.6),
             ("LINEBELOW", (0, 0), (-1, -1), 0.3, RULE), ("BOX", (0, 0), (-1, -1), 0.5, RULE)]
    for i in range(1, len(data)):
        if i % 2 == 0:
            style.append(("BACKGROUND", (0, i), (-1, i), colors.HexColor("#f7f6f2")))
    t.setStyle(TableStyle(style))
    out = [t]
    if caption:
        out.append(Paragraph(M(caption), ST["cap"]))
    else:
        out.append(Spacer(1, 2.5 * mm))
    return [KeepTogether(out)] if keep else out


def CARD(sid, title, goal, decisions, depends=None):
    """decisions: list of dicts with keys d (decision), why, risk (optional), check, fail (optional)."""
    lab_w = 19 * mm
    cw = [lab_w, WIDTH - lab_w]
    base = [("VALIGN", (0, 0), (-1, -1), "TOP"), ("LEFTPADDING", (0, 0), (-1, -1), 4), ("RIGHTPADDING", (0, 0), (-1, -1), 4),
            ("TOPPADDING", (0, 0), (-1, -1), 2.4), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.8),
            ("LINEBEFORE", (0, 0), (0, -1), 0.6, NAVY), ("LINEAFTER", (-1, 0), (-1, -1), 0.6, NAVY)]
    head = Table([[Paragraph(M(f"{sid}   {title}"), ST["cardhead"]), ""],
                  [Paragraph("Goal", ST["cardlab"]), Paragraph(M(goal), ST["card"])]], colWidths=cw)
    head.setStyle(TableStyle(base + [("SPAN", (0, 0), (1, 0)), ("BACKGROUND", (0, 0), (-1, 0), NAVY)]))
    blocks = []                                  # one table per decision, so that a decision is never split over pages
    for k, dct in enumerate(decisions, 1):
        data = [[Paragraph(f"Decision {k}", ST["cardlab"]), Paragraph("<b>" + M(dct["d"]) + "</b>", ST["card"])]]
        style = base + [("BACKGROUND", (0, 0), (-1, 0), PALE), ("LINEABOVE", (0, 0), (-1, 0), 0.6, NAVY)]
        for key, lab in (("why", "Why"), ("risk", "Risk"), ("check", "Check"), ("fail", "If it fails")):
            if not dct.get(key):
                continue
            data.append([Paragraph(lab, ST["cardlab"]), Paragraph(M(dct[key]), ST["card"])])
            if key == "check":
                style.append(("BACKGROUND", (0, len(data) - 1), (-1, len(data) - 1), WARM))
        if k == len(decisions):
            if depends:
                data.append([Paragraph("Depends on", ST["cardlab"]), Paragraph(M(depends), ST["card"])])
                style.append(("LINEABOVE", (0, len(data) - 1), (-1, len(data) - 1), 0.6, NAVY))
            style.append(("LINEBELOW", (0, len(data) - 1), (-1, len(data) - 1), 0.6, NAVY))
        t = Table(data, colWidths=cw, splitByRow=1)
        t.setStyle(TableStyle(style))
        blocks.append(t)
    return ([CondPageBreak(45 * mm), KeepTogether([head, blocks[0]])] + [KeepTogether([b]) for b in blocks[1:]]
            + [Spacer(1, 4 * mm)])


def BOX(title, items, numbered=True):
    rows = [[Paragraph(M(title), ParagraphStyle("bt", parent=ST["card"], fontName="A-B", fontSize=9.2, textColor=NAVY))]]
    for i, it in enumerate(items, 1):
        txt = (f"**{i}.** " if numbered else "") + it
        rows.append([Paragraph(M(txt), ParagraphStyle("bi", parent=ST["card"], fontSize=8.6, leading=11.2))])
    t = Table(rows, colWidths=[WIDTH])
    t.setStyle(TableStyle([("BACKGROUND", (0, 0), (-1, -1), PALE), ("BOX", (0, 0), (-1, -1), 0.6, BLUE),
                           ("LEFTPADDING", (0, 0), (-1, -1), 6), ("RIGHTPADDING", (0, 0), (-1, -1), 6),
                           ("TOPPADDING", (0, 0), (-1, -1), 2.2), ("BOTTOMPADDING", (0, 0), (-1, -1), 2.4)]))
    return [t, Spacer(1, 3.5 * mm)]


def FIG(path, caption, width_frac=1.0):
    from PIL import Image as PI
    with PI.open(path) as im:
        w, h = im.size
    W = WIDTH * width_frac
    return [KeepTogether([Image(path, width=W, height=W * h / w), Paragraph(M(caption), ST["cap"])])]


def HFIG(level, title, path, caption, width_frac=1.0, intro=None):
    """Heading, optional intro paragraph and figure kept on one page, so a heading never stands alone."""
    from PIL import Image as PI
    with PI.open(path) as im:
        w, h = im.size
    W = WIDTH * width_frac
    parts = [Paragraph(M(title), ST["h1" if level == 1 else "h2"])]
    if intro:
        parts.append(Paragraph(M(intro), ST["body"]))
    parts += [Image(path, width=W, height=W * h / w), Paragraph(M(caption), ST["cap"])]
    return [KeepTogether(parts)]


def _page(canvas, doc):
    canvas.saveState()
    canvas.setFont("A", 7.3); canvas.setFillColor(MUTED)
    canvas.drawString(LM, 9.5 * mm, doc.footer_text)
    canvas.drawRightString(PAGE_W - RM, 9.5 * mm, f"Page {doc.page}")
    canvas.setStrokeColor(RULE); canvas.setLineWidth(0.4)
    canvas.line(LM, 12.5 * mm, PAGE_W - RM, 12.5 * mm)
    canvas.restoreState()


def check_glyphs():
    """Every character used must exist in the body font, otherwise it would print as a box."""
    from fontTools.ttLib import TTFont as FT
    cmap = FT(f"{FONTS}/arial.ttf").getBestCmap()
    missing = sorted({ch for t in USED_TEXT for ch in t if ord(ch) > 126 and ord(ch) not in cmap})
    return missing


def build(path, story, title, footer):
    flat = []
    for s in story:
        flat.extend(s if isinstance(s, list) else [s])
    doc = SimpleDocTemplate(path, pagesize=A4, leftMargin=LM, rightMargin=RM, topMargin=TM, bottomMargin=BM,
                            title=title, author="Ogallala groundwater buffering project", subject="Detailed analysis plan")
    doc.footer_text = footer
    doc.build(flat, onFirstPage=_page, onLaterPages=_page)
    return doc.page
