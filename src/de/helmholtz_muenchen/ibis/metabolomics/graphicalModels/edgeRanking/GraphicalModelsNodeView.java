package de.helmholtz_muenchen.ibis.metabolomics.graphicalModels.edgeRanking;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeView</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class GraphicalModelsNodeView extends RNodeView<GraphicalModelsNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link GraphicalModelsNodeModel})
     */
    protected GraphicalModelsNodeView(final GraphicalModelsNodeModel nodeModel) {
        super(nodeModel);
    }
}

