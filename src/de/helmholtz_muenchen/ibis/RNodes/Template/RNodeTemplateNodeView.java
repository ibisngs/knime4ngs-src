package de.helmholtz_muenchen.ibis.RNodes.Template;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeView</code> for the "RNodeTemplate" Node.
 * This is a template for a node, which executes a R-Script.
 *
 * @author Jonas Zierer
 */
public class RNodeTemplateNodeView extends RNodeView<RNodeTemplateNodeModel> {

    /**
     * Creates a new view.
     * 
     * @param nodeModel The model (class: {@link RNodeTemplateNodeModel})
     */
    protected RNodeTemplateNodeView(final RNodeTemplateNodeModel nodeModel) {
        super(nodeModel);
    }
}

