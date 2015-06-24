package de.helmholtz_muenchen.ibis.ngs.grouplofgenes;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GroupLoFGenes" Node.
 * 
 *
 * @author tim.jeske
 */
public class GroupLoFGenesNodeFactory 
        extends NodeFactory<GroupLoFGenesNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GroupLoFGenesNodeModel createNodeModel() {
        return new GroupLoFGenesNodeModel();
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
    public NodeView<GroupLoFGenesNodeModel> createNodeView(final int viewIndex,
            final GroupLoFGenesNodeModel nodeModel) {
        return new GroupLoFGenesNodeView(nodeModel);
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
        return new GroupLoFGenesNodeDialog();
    }

}

