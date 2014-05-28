package de.helmholtz_muenchen.ibis.plotting.barplot;

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
public class BarplotNodeFactory extends NodeFactory<BarplotNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public BarplotNodeModel createNodeModel() {
        return new BarplotNodeModel();
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
    public NodeView<BarplotNodeModel> createNodeView(final int viewIndex, final BarplotNodeModel nodeModel) {
        return new RNodeView<BarplotNodeModel>(nodeModel);
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
        return new BarplotNodeDialog();
    }

}

