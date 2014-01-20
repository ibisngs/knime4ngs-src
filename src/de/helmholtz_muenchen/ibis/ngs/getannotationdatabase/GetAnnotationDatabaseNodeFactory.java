package de.helmholtz_muenchen.ibis.ngs.getannotationdatabase;

import org.knime.core.node.NodeDialogPane;
import org.knime.core.node.NodeFactory;
import org.knime.core.node.NodeView;

/**
 * <code>NodeFactory</code> for the "GetAnnotationDatabase" Node.
 * 
 *
 * @author 
 */
public class GetAnnotationDatabaseNodeFactory 
        extends NodeFactory<GetAnnotationDatabaseNodeModel> {

    /**
     * {@inheritDoc}
     */
    @Override
    public GetAnnotationDatabaseNodeModel createNodeModel() {
        return new GetAnnotationDatabaseNodeModel();
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
    public NodeView<GetAnnotationDatabaseNodeModel> createNodeView(final int viewIndex,
            final GetAnnotationDatabaseNodeModel nodeModel) {
        return new GetAnnotationDatabaseNodeView(nodeModel);
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
        return new GetAnnotationDatabaseNodeDialog();
    }

}

