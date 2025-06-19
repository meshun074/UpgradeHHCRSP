package org.example.GA;

import org.example.Data.InstancesClass;
import org.example.Main;

import java.util.ArrayList;
import java.util.Collections;
import java.util.List;

public class Population {
    public static List<Chromosome> initialize(int popSize, int chLength) {
        List<Integer> patients= new ArrayList<>();
        List<Integer> patientOrder;
        for (int s = 0; s < chLength; s++) {
            patients.add(s);
        }
        InstancesClass instances = Main.instance;
        List<RouteInitializer> tempPopulation = new ArrayList<>(popSize);
        AssignPatients as = new AssignPatients(instances);
        int caregiverNumber = instances.getCaregivers().length;
        for (int i = 0; i < popSize; i++) {
            patientOrder = new ArrayList<>(patients);
            Collections.shuffle(patientOrder);
            RouteInitializer ri = new RouteInitializer(patientOrder,caregiverNumber,0.0);
            tempPopulation.add(ri);
        }
        as.Evaluate(tempPopulation);

        List<Chromosome> population = new ArrayList<>(popSize);
        for (RouteInitializer ri : tempPopulation) {
            Chromosome ch = new Chromosome(ri.getCaregiversRoute(),ri.getSolutionCost(),false, ri.getCaregiversShift());
            population.add(ch);
            population.getLast().setTotalTravelCost(ri.getTotalTravelCost());
            population.getLast().setTotalTardiness(ri.getTotalTardiness());
            population.getLast().setHighestTardiness(ri.getHighestTardiness());
        }
        return population;
    }
}
