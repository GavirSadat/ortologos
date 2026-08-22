from ete3 import Tree, TreeStyle, SeqMotifFace, TextFace, NodeStyle
import os

tree_file = "alineamiento_limpio.fasta.treefile"
alignment_file = "alineamiento_limpio.fasta"
coords_file = "coordenadas_della.txt"

if os.path.exists(tree_file) and os.path.exists(alignment_file) and os.path.exists(coords_file):
    # 1. Cargar topología
    tree = Tree(tree_file, format=1)
    
    # ---> NOVEDAD: Enraizar el árbol usando la licofita como grupo externo basal
    tree.set_outgroup("Selaginella_moellendorffii")

    # 2. Cargar MSA curado para mapear las coordenadas homólogas
    seqs = {}
    with open(alignment_file, "r") as f:
        name = ""
        for line in f:
            if line.startswith(">"):
                name = line.strip()[1:] # IQ-TREE trunca nombres largos
                seqs[name] = ""
            else:
                seqs[name] += line.strip()

    # 3. Recuperar coordenadas homólogas exactas del Módulo 5
    with open(coords_file, "r") as f:
        start_col, end_col = map(int, f.read().split(","))

    motifs = {}
    for name, seq in seqs.items():
        motifs[name] = [[start_col, end_col, "[]", None, 10, "red", "red", "arial|8|white|Región_DELLA"]]

    # 4. Integrar Secuencia-Alineamiento-Topología
    for leaf in tree.iter_leaves():
        if leaf.name in seqs:
            leaf.add_face(SeqMotifFace(seqs[leaf.name], motifs=motifs.get(leaf.name, []), seq_format="()"), column=1, position="aligned")

    # 5. Integración de Soportes Estadísticos (UFBoot) con colores forzados
    for node in tree.traverse():
        if not node.is_leaf():
            support_val = 0
            if node.name and "/" in node.name:
                try:
                    support_val = int(float(node.name.split("/")[1]))
                except ValueError:
                    support_val = 0
            elif hasattr(node, 'support') and node.support:
                support_val = int(node.support)
                
            nstyle = NodeStyle()
            nstyle["shape"] = "sphere" # Forzar el dibujo del nodo
            nstyle["size"] = 9         # Aumentar tamaño para mejor visibilidad
            
            # Esquema semafórico
            if support_val >= 90:
                nstyle["fgcolor"] = "green"
            elif support_val >= 70:
                nstyle["fgcolor"] = "orange"
            else:
                nstyle["fgcolor"] = "red"
            
            node.set_style(nstyle)
            # Anotación del valor de soporte
            node.add_face(TextFace(f" {support_val} ", fsize=8, fgcolor="black"), column=0, position="branch-top")

    # 6. Renderizado de la figura integrada
    ts = TreeStyle()
    ts.show_leaf_name = True
    ts.draw_guiding_lines = True
    ts.title.add_face(TextFace("Análisis Filogenómico Integrado: Evolución del Dominio DELLA \n", fsize=14, bold=True), column=0)
    
    # 7. Guardar en la carpeta results
    pdf_path = "../results/Integracion_DELLA_ETE3.pdf"
    os.makedirs(os.path.dirname(pdf_path), exist_ok=True)
    
    tree.render(pdf_path, tree_style=ts, w=1200)
    print(f"[OK] Diagrama de integración (ENRAIZADO) guardado en: {pdf_path}")
else:
    print("Error: Archivos requeridos no encontrados.")
