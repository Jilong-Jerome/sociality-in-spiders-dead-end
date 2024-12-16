#!/bin/env python3
import sys, os
sys.path.insert(0, os.path.realpath('/faststorage/project/EcoGenetics/people/Jeppe_Bayer/population_genetics/scripts/mapping/workflow_source'))	
from workflow_source import *

gwf = mapping_resequencing_data_population_genetics_workflow()
