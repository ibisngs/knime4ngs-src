package de.helmholtz_muenchen.icb.epigenreg.macs2;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "Macs2" Node.
 * 
 *
 * @author 
 */
public class Macs2NodeFactory 
        extends NodeFactory<Macs2NodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public Macs2NodeModel createNodeModel() {
        return new Macs2NodeModel();
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
    public NodeView<Macs2NodeModel> createNodeView(final int viewIndex,
            final Macs2NodeModel nodeModel) {
        return new Macs2NodeView(nodeModel);
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
        return new Macs2NodeDialog();
    }

}

