package de.helmholtz_muenchen.ibis.ngs.gatkexcludevariants;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GATKExcludeVariants" Node.
 * 
 *
 * @author Tim Jeske
 */
public class GATKExcludeVariantsNodeFactory 
        extends NodeFactory<GATKExcludeVariantsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GATKExcludeVariantsNodeModel createNodeModel() {
        return new GATKExcludeVariantsNodeModel();
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
    public NodeView<GATKExcludeVariantsNodeModel> createNodeView(final int viewIndex,
            final GATKExcludeVariantsNodeModel nodeModel) {
        return new GATKExcludeVariantsNodeView(nodeModel);
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
        return new GATKExcludeVariantsNodeDialog();
    }

}

