package de.helmholtz_muenchen.ibis.plotting.scatterplot;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

import de.helmholtz_muenchen.ibis.utils.abstractNodes.RNode.RNodeView;

/**
 * <code>NodeFactory</code> for the "Boxplot" Node.
 * 
 *
 * @author Jonas Zierer
 */
public class ScatterplotNodeFactory extends NodeFactory<ScatterplotNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public ScatterplotNodeModel createNodeModel() {
        return new ScatterplotNodeModel();
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public int getNrNodeViews() {
        return 1;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeView<ScatterplotNodeModel> createNodeView(final int viewIndex, final ScatterplotNodeModel nodeModel) {
        return new RNodeView<ScatterplotNodeModel>(nodeModel);
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public boolean hasDialog() {
        return true;
    }

    /**
     * {@inheritDoc}
     */
    @Override
    public NodeDialogPane createNodeDialogPane() {
        return new ScatterplotNodeDialog();
    }

}

