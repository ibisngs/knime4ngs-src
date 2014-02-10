package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.modelExtraction;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeView</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class GraphicalModelExtractionNodeView extends RNodeView<GraphicalModelExtractionNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link GraphicalModelExtractionNodeModel})
     */
    protected GraphicalModelExtractionNodeView(final GraphicalModelExtractionNodeModel nodeModel) {
        super(nodeModel);
    }
}

