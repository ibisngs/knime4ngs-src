package de.helmholtz_muenchen.ibis.ngs.filtergatkvariants;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "FilterGATKVariants" Node.
 * 
 *
 * @author Maximilian Hastreiter
 */
public class FilterGATKVariantsNodeFactory 
        extends NodeFactory<FilterGATKVariantsNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public FilterGATKVariantsNodeModel createNodeModel() {
        return new FilterGATKVariantsNodeModel();
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
    public NodeView<FilterGATKVariantsNodeModel> createNodeView(final int viewIndex,
            final FilterGATKVariantsNodeModel nodeModel) {
        return new FilterGATKVariantsNodeView(nodeModel);
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
        return new FilterGATKVariantsNodeDialog();
    }

}

