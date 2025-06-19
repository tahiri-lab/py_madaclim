from graphviz import Digraph

# Create the main data flow diagram using Graphviz
dot = Digraph(comment='Coffea Caffeine Prediction Pipeline', format='png')

# Main pipeline stages
dot.node('A', 'Raw Data Sources\n(CIRAD coordinates,\nGenomic species,\nMadaclim rasters)')
dot.node('B', 'Data Integration\n+ Correction\n- Geospatial merge\n- Environmental extraction\n- Jan anomaly fix')
dot.node('C', 'Preprocessing\n- Categorical encoding\n- Imputation\n- Filtering (valid coords)')
dot.node('D', 'Feature Selection\n- Correlation Clustering\n- Variance / RFC criteria')
dot.node('E', 'Stage 1: Classification\n(Random Forest Classifier)\nPredict caffeine presence')
dot.node('F', 'Stage 2: Regression\n(Random Forest Regressor)\nPredict caffeine content')
dot.node('G', 'Model Evaluation\n- Accuracy, AUC, F1\n- MSE, R²')
dot.node('H', 'Visualization\n- Phylogenetic Map\n- Caffeine Gradient Mapping\n- Feature Importance')
dot.node('I', 'Final Outputs\n- CSV + plots\n- Feature table\n- Integrated maps')

# Connections
dot.edges(['AB', 'BC', 'CD', 'DE', 'EF', 'FG', 'GH', 'HI'])

# Display the diagram
dot.render('../images/coffea_pipeline', view=False)

