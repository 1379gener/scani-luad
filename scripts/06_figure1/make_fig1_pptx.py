# Figure 1 可编辑版（PPTX）：与 make_fig1.py 同一布局、同一文字；坐标系 100 x 72 → 7.2 x 5.0 英寸
from pptx import Presentation
from pptx.util import Inches, Pt
from pptx.dml.color import RGBColor
from pptx.enum.shapes import MSO_SHAPE, MSO_CONNECTOR
from pptx.enum.text import PP_ALIGN, MSO_ANCHOR
from pptx.oxml.ns import qn
from lxml import etree
W, H = 7.2, 5.0
prs = Presentation(); prs.slide_width = Inches(W); prs.slide_height = Inches(H)
s = prs.slides.add_slide(prs.slide_layouts[6])
X = lambda x: Inches(x / 100 * W); Y = lambda y: Inches((72 - y) / 72 * H)   # matplotlib y 向上 → pptx y 向下
INK, MUTED = RGBColor(0x22, 0x22, 0x22), RGBColor(0x66, 0x66, 0x66)
FILL_MAIN, FILL_OUT, FILL_SIDE = RGBColor(0xEE, 0xF3, 0xF8), RGBColor(0xF6, 0xF6, 0xF6), RGBColor(0xFF, 0xFF, 0xFF)
EDGE_MAIN, EDGE_OUT = RGBColor(0x4A, 0x6A, 0x8A), RGBColor(0x7A, 0x7A, 0x7A)
def box(x, y, w, h, title, body, fill, edge, dashed=False):
    sh = s.shapes.add_shape(MSO_SHAPE.ROUNDED_RECTANGLE, X(x), Y(y + h), X(w) - X(0), Y(y) - Y(y + h))
    sh.adjustments[0] = 0.08
    sh.fill.solid(); sh.fill.fore_color.rgb = fill
    sh.line.color.rgb = edge; sh.line.width = Pt(1.1)
    if dashed:
        ln = sh.line._get_or_add_ln(); d = etree.SubElement(ln, qn("a:prstDash")); d.set("val", "dash")
    tf = sh.text_frame; tf.word_wrap = True; tf.vertical_anchor = MSO_ANCHOR.TOP
    for m in ("margin_left", "margin_right"): setattr(tf, m, Inches(0.04))
    tf.margin_top = Inches(0.06)
    p = tf.paragraphs[0]; p.alignment = PP_ALIGN.CENTER
    r = p.add_run(); r.text = title; r.font.size = Pt(8.6); r.font.bold = True; r.font.color.rgb = INK; r.font.name = "Arial"
    for line in body.split("\n"):
        q = tf.add_paragraph(); q.alignment = PP_ALIGN.CENTER
        r = q.add_run(); r.text = line; r.font.size = Pt(7.2); r.font.color.rgb = MUTED; r.font.name = "Arial"
def arrow(x1, y1, x2, y2):
    c = s.shapes.add_connector(MSO_CONNECTOR.STRAIGHT, X(x1), Y(y1), X(x2), Y(y2))
    c.line.color.rgb = RGBColor(0x55, 0x55, 0x55); c.line.width = Pt(1.0)
    ln = c.line._get_or_add_ln(); t = etree.SubElement(ln, qn("a:tailEnd")); t.set("type", "triangle")
def note(x, y, text):
    tb = s.shapes.add_textbox(X(x - 9), Y(y + 4), X(18) - X(0), Y(y - 4) - Y(y + 4)); tf = tb.text_frame; tf.word_wrap = True
    for i, line in enumerate(text.split("\n")):
        p = tf.paragraphs[0] if i == 0 else tf.add_paragraph(); p.alignment = PP_ALIGN.CENTER
        r = p.add_run(); r.text = line; r.font.size = Pt(6.8); r.font.color.rgb = MUTED; r.font.name = "Arial"
box(2, 52, 44, 17, "Single-cell RNA-seq (GSE131907)",
    "Fibroblasts as annotated by the original study\n(primary tumour, normal lung, metastases)", FILL_MAIN, EDGE_MAIN)
box(54, 52, 44, 17, "myCAF signature (50 genes)",
    "Myofibroblasts vs. other fibroblast states\n(FindAllMarkers; top 50 by log2 fold change)", FILL_MAIN, EDGE_MAIN)
arrow(46.5, 60.5, 53.5, 60.5)
box(27, 30, 50, 15, "SCANI", "Z(myCAF ssGSEA) + Z(TIMP1) − Z(PDPN)\neach component standardized within cohort", FILL_MAIN, EDGE_MAIN)
arrow(76, 51.5, 62, 45.6)
note(86, 47.5, "TIMP1 and PDPN\nselected in\nTCGA-LUAD")
box(2, 3, 30, 19, "TCGA-LUAD", "Discovery, stage II–IV\nCox (per SD), adjusted Cox,\ntime-dependent ROC, KM", FILL_OUT, EDGE_OUT)
box(35, 3, 30, 19, "GSE72094", "Independent validation\n(microarray platform)\nCox (per SD), adjusted Cox, KM", FILL_OUT, EDGE_OUT)
box(68, 3, 30, 19, "CPTAC-LUAD", "Association with stage\n(Spearman correlation\nwith pathological stage)", FILL_OUT, EDGE_OUT)
for xc in (17, 50, 83): arrow(50, 29.5, xc, 22.6)
box(2, 30, 21, 15, "Visium (GSE307534)", "Descriptive stromal\nprogram; not a\nSCANI component", FILL_SIDE, EDGE_OUT, dashed=True)
prs.save(__import__("os").path.join(__import__("os").environ.get("FIG_OUT", "figures"), "Figure1_editable.pptx")); print("saved")
