package org.example.GA;

import org.example.Data.InstancesClass;
import org.example.Data.Patient;
import org.example.Data.Required_Caregiver;

import java.util.*;

import static org.example.GA.EvaluationFunction.removeAffectedPatient;
import static org.example.GA.GeneticAlgorithm.conflictCheck;

public class BestCostRouteCrossoverNewNew implements Runnable {
    @Override
    public void run() {
        ga.getCrossoverChromosomes().add(Crossover());
    }

    private final GeneticAlgorithm ga;
    private final int r;
    private final int index;
    private final Chromosome p1, p2;
    private static Patient[] allPatients;
    private static InstancesClass dataset;
    private static int allCaregivers;
    private static double[][] distances;
    private final Random rand;
    private final int[] routeEndPoint;
    private final double[] routesCurrentTime;
    private final double[] highestAndTotalTardiness;
    private final Set<Integer> track;

    public BestCostRouteCrossoverNewNew(GeneticAlgorithm ga, int r, Chromosome p1, Chromosome p2, int index, Random rand) {
        this.ga = ga;
        this.r = r;
        this.index = index;
        this.p1 = p1;
        this.p2 = p2;
        this.rand = rand;
        this.routeEndPoint = new int[allCaregivers];
        this.routesCurrentTime = new double[allCaregivers];
        this.highestAndTotalTardiness = new double[2];
        this.track = new HashSet<>(100);
    }

    public static void initialize(InstancesClass instance) {
        allPatients = instance.getPatients();
        distances = instance.getDistances();
        allCaregivers = instance.getCaregivers().length;
        dataset = instance;
    }

    public Chromosome Crossover() {
        // Initialize variables
        Set<Integer> selectRoute = new HashSet<>(p2.getGenes()[r]);
        List<Integer>[] p1Routes = p1.getGenes();
        List<Integer>[] c1Routes = new ArrayList[p1Routes.length];

        //Remove patients of selected route from parents
        for (int i = 0; i < p1Routes.length; i++) {
            List<Integer> route = new ArrayList<>(p1Routes[i].size());
            for (int j = 0; j < p1Routes[i].size(); j++) {
                int patient = p1Routes[i].get(j);
                if (!selectRoute.contains(patient)) {
                    route.add(patient);
                }
            }
            c1Routes[i] = route;
        }
        List<Integer> route = new ArrayList<>(selectRoute);
//        Collections.shuffle(route);
        Chromosome cTemp = new Chromosome(c1Routes, 0.0, false);
        EvaluationFunction.Evaluate(cTemp);
//        System.out.println("removed: "+route);
//        System.out.println("cTempo: "+cTemp+" - "+cTemp.getFitness());
        for (int i = 0; i < route.size(); i++) {
            int patient = route.get(i);
            Patient p = allPatients[patient];
            double bestCost = Double.MAX_VALUE;
            int bestFirst = -1;
            int bestSecond = -1;
            int bestM = -1;
            int bestN = -1;
            cTemp.buildPatientRouteMap();
            Shift[] shifts = cTemp.getCaregiversRouteUp();
            boolean isInvalid = cTemp.getFitness() == Double.POSITIVE_INFINITY;
            List<CaregiverPair> caregiverPairs = p.getAllPossibleCaregiverCombinations();
            if (p.getRequired_caregivers().length > 1) {
                for (int x = 0; x < caregiverPairs.size(); x++) {
                    CaregiverPair caregiverPair = caregiverPairs.get(x);
                    int first = caregiverPair.getFirst();
                    int second = caregiverPair.getSecond();
//                    System.out.println("first: "+first+" second: "+second);
                    for (int m = 0; m <= c1Routes[first].size(); m++) {
                        for (int n = 0; n <= c1Routes[second].size(); n++) {
//                            System.out.println("cTemp: m "+m+" - n"+n+" - "+cTemp);
                            if (noEvaluationConflicts(c1Routes[first], c1Routes[second], m, n)) {
                                double tempCost = calMoveCost(first, m, second, n, patient, cTemp, bestCost, shifts, isInvalid);
//                                System.out.println("tempCost: "+tempCost);
                                if (bestCost == Double.MAX_VALUE || bestCost - tempCost > 0.001 && bestCost != Double.POSITIVE_INFINITY && tempCost != Double.POSITIVE_INFINITY || tempCost <= bestCost && rand.nextBoolean()) {
                                    bestCost = tempCost;
                                    bestFirst = first;
                                    bestSecond = second;
                                    bestM = m;
                                    bestN = n;
                                }
                            }
                        }
                    }
                }
                if (bestCost != Double.MAX_VALUE) {
                    c1Routes = cTemp.getGenes();
                    c1Routes[bestFirst].add(bestM, patient);
                    c1Routes[bestSecond].add(bestN, patient);
                    EvaluationFunction.Evaluate(cTemp);
                } else {
                    System.out.println("no route found");
                }
            } else {
                for (int x = 0; x < caregiverPairs.size(); x++) {
                    CaregiverPair caregiverPair = caregiverPairs.get(x);
                    int first = caregiverPair.getFirst();
                    for (int k = 0; k <= c1Routes[first].size(); k++) {
                        double tempCost = calMoveCost(first,k,-1,-1,patient,cTemp,bestCost,shifts,isInvalid);
                        if (bestCost == Double.MAX_VALUE || bestCost - tempCost > 0.001 && bestCost != Double.POSITIVE_INFINITY && tempCost != Double.POSITIVE_INFINITY || tempCost <= bestCost && rand.nextBoolean()) {
                            bestCost = tempCost;
                            bestFirst = first;
                            bestM = k;
                        }
                    }
                }
                if (bestCost != Double.MAX_VALUE) {
                    c1Routes = cTemp.getGenes();
                    c1Routes[bestFirst].add(bestM, patient);
                    EvaluationFunction.Evaluate(cTemp);
                } else {
                    System.out.println("no route found " + caregiverPairs.size());
                }
            }
//            System.out.println("cTempo: "+cTemp+" - "+cTemp.getFitness());
        }
        if (Math.round(cTemp.getFitness()) == Math.round(157.33566666666667)) {
            System.exit(1);
        }
        cTemp.setCrossIndex(index);
        return cTemp;
    }

    private double calMoveCost(int first, int m, int second, int n, int patient, Chromosome cTemp, double bestCost, Shift[] shifts, boolean isInvalid) {
        Arrays.fill(routeEndPoint, -1);
        if (isInvalid) {
            List<Integer>[] c1Routes = cTemp.getGenes();
            c1Routes[first].add(m, patient);
            if(second !=-1){
                c1Routes[second].add(n, patient);
            }
            Chromosome temp = new Chromosome(c1Routes, 0.0, false);
            c1Routes[first].remove(Integer.valueOf(patient));
            if(second !=-1){
                c1Routes[second].remove(Integer.valueOf(patient));
            }
            EvaluationFunction.EvaluateNew(temp);
            return temp.getFitness();
        }
        int size = 1;
        routeEndPoint[first] = m;
        if (second != -1) {
            size++;
            routeEndPoint[second] = n;
        }

        int[] routeMove = new int[size];
        int[] positionMove = new int[size];
        routeMove[0] = first;
        positionMove[0] = m;
        if (size > 1) {
            routeMove[1] = second;
            positionMove[1] = n;
        }
        removeAffectedPatient(routeMove, positionMove, cTemp, routeEndPoint);

        double totalTravelCost = 0.0;
        highestAndTotalTardiness[0] = 0;
        highestAndTotalTardiness[1] = 0;
        for (int i = 0; i < routeEndPoint.length; i++) {
            List<Double> currentTime = shifts[i].getCurrentTime();
            List<Double> travelCost = shifts[i].getTravelCost();
            List<Double> tardiness = shifts[i].getTardiness();
            List<Double> maxTardiness = new ArrayList<>(shifts[i].getMaxTardiness());

            if (routeEndPoint[i] != -1) {
                int index = routeEndPoint[i];
                routesCurrentTime[i] = currentTime.get(index);
                highestAndTotalTardiness[0] = Math.max(maxTardiness.get(index), highestAndTotalTardiness[0]);
                highestAndTotalTardiness[1] += tardiness.get(index);
            } else {
                highestAndTotalTardiness[0] = Math.max(maxTardiness.get(maxTardiness.size() - 1), highestAndTotalTardiness[0]);
                highestAndTotalTardiness[1] += tardiness.get(tardiness.size() - 1);
            }
            if (i == first || i == second) {
                int index = routeEndPoint[i];
                totalTravelCost += travelCost.get(index);
            } else {
                totalTravelCost += travelCost.get(travelCost.size() - 1);
            }
        }


        //Distance calculation
        List<Integer>[] genes = cTemp.getGenes();
        genes[first].add(m, patient);
        if(second !=-1){
            genes[second].add(n, patient);
        }
        for (int i = 0; i < routeEndPoint.length; i++) {
            List<Integer> route = genes[i];
            int routeEnd = routeEndPoint[i];
            if (i == first || i == second) {
                for (int j = routeEnd; j <= route.size(); j++) {
                    if (j == 0) {
                        int nextIndex = route.get(j)+1;
                        totalTravelCost += distances[0][nextIndex];
                    } else if (j == route.size()) {
                        int prevIndex = route.get(j-1)+1;
                        totalTravelCost += distances[prevIndex][0];
                    } else {
                        int nextIndex = route.get(j)+1;
                        int prevIndex = route.get(j-1)+1;
                        totalTravelCost += distances[prevIndex][nextIndex];
                    }

                }
            }
        }
        totalTravelCost = 1 / 3d * totalTravelCost;

        double solutionCost = totalTravelCost + (1 / 3d * highestAndTotalTardiness[0]) + (1 / 3d * highestAndTotalTardiness[1]);
        if (solutionCost > bestCost) {
            genes[first].remove(Integer.valueOf(patient));
            if(second !=-1){
                genes[second].remove(Integer.valueOf(patient));
            }
            return solutionCost;
        }

        track.clear();
        //Tardiness calculation
        for (int i = 0; i < routeEndPoint.length; i++) {
            List<Integer> route = genes[i];
            int routeEnd = routeEndPoint[i];
            if (routeEnd != -1) {
                for (int j = routeEnd; j < route.size(); j++) {
                    if (j == 0) {
                        solutionCost = patientIsAssigned(genes, i, -1, route.get(j), totalTravelCost, routesCurrentTime, highestAndTotalTardiness, track);
                    } else {
                        solutionCost = patientIsAssigned(genes, i,  route.get(j - 1), route.get(j), totalTravelCost, routesCurrentTime, highestAndTotalTardiness, track);
                    }
                    if (solutionCost == Double.POSITIVE_INFINITY || solutionCost > bestCost) {
                        genes[first].remove(Integer.valueOf(patient));
                        if(second !=-1){
                            genes[second].remove(Integer.valueOf(patient));
                        }
                        return solutionCost;
                    }

                    track.clear();
                }
            }
        }
        genes[first].remove(Integer.valueOf(patient));
        if(second !=-1){
            genes[second].remove(Integer.valueOf(patient));
        }

        return solutionCost;
    }

    private double patientIsAssigned(List<Integer>[] genes, int route1, int curPatient1, int patient, double totalTravelCost, double[] routesCurrentTime, double[] highestAndTotalTardiness, Set<Integer> track) {
        Patient p = allPatients[patient];
        double[] timeWindow = p.getTime_window();
        int currentLocation1 = curPatient1 + 1;
        int nextLocation = patient + 1;

        double arrivalTime1 = routesCurrentTime[route1] + distances[currentLocation1][nextLocation];
        double startTime1 = Math.max(arrivalTime1, timeWindow[0]);

        if (p.getRequired_caregivers().length > 1) {
            if (!track.add(patient)) {
                return Double.POSITIVE_INFINITY;
            }

            int route2 = findOtherCaregiver(patient, route1, genes, totalTravelCost,  routesCurrentTime, highestAndTotalTardiness, track);
            if (route2 > allCaregivers - 1) {
                return Double.POSITIVE_INFINITY;
            }

            Required_Caregiver[] requiredCaregivers = p.getRequired_caregivers();

            int position2 = routeEndPoint[route2] - 1;
            int curPatient2;
            if (position2 == -1) {
                curPatient2 = -1;
            }else {
                curPatient2 = genes[route2].get(position2);
            }

            if (SwapRoutes(p, route1, route2)) {
                int temp = route1;
                route1 = route2;
                route2 = temp;
                temp = curPatient1;
                curPatient1 = curPatient2;
                curPatient2 = temp;
                currentLocation1 = curPatient1 + 1;
                arrivalTime1 = routesCurrentTime[route1] + distances[currentLocation1][nextLocation];
                startTime1 = Math.max(arrivalTime1, timeWindow[0]);
            }

            int currentLocation2 = curPatient2 + 1;
            double arrivalTime2 = routesCurrentTime[route2] + distances[currentLocation2][nextLocation];
            double startTime2 = Math.max(arrivalTime2, timeWindow[0]);

            if (p.getSynchronization().getType().equals("sequential")) {
                double[] syncDistances = p.getSynchronization().getDistance();
                startTime2 = Math.max(startTime2, startTime1 + syncDistances[0]);
                if (startTime2 - startTime1 > syncDistances[1]) {
                    startTime1 = startTime2 - syncDistances[1];
                }
                double tardiness = Math.max(0, startTime1 - timeWindow[1]);
                highestAndTotalTardiness[0] = Math.max(highestAndTotalTardiness[0], tardiness);
                highestAndTotalTardiness[1] += tardiness;
                tardiness = Math.max(0, startTime2 - timeWindow[1]);
                highestAndTotalTardiness[0] = Math.max(highestAndTotalTardiness[0], tardiness);
                highestAndTotalTardiness[1] += tardiness;
            }else {
                double startTime = Math.max(startTime1, startTime2);
                double tardiness = Math.max(0, startTime - timeWindow[1]);
                highestAndTotalTardiness[0] = Math.max(highestAndTotalTardiness[0], tardiness);
                highestAndTotalTardiness[1] += (2*tardiness);
                startTime1 = startTime2 = startTime;
            }
            routeEndPoint[route1]++;
            routeEndPoint[route2]++;
            routesCurrentTime[route1] = startTime1 + p.getRequired_caregivers()[0].getDuration();
            routesCurrentTime[route2] = startTime2 + p.getRequired_caregivers()[1].getDuration();

        } else {
            double tardiness = Math.max(0, startTime1 - timeWindow[1]);
            highestAndTotalTardiness[0] = Math.max(highestAndTotalTardiness[0], tardiness);
            highestAndTotalTardiness[1] += tardiness;
            routeEndPoint[route1]++;
            routesCurrentTime[route1] = startTime1 + p.getRequired_caregivers()[0].getDuration();
        }
        return totalTravelCost + (1 / 3d * highestAndTotalTardiness[0]) + (1 / 3d * highestAndTotalTardiness[1]);
    }

    private int findOtherCaregiver(int patient, int route1, List<Integer>[] genes, double totalTravelCost, double[] routesCurrentTime, double[] highestAndTotalTardiness, Set<Integer> track) {
        for (int i = 0; i < genes.length; i++) {
            if (i != route1 && genes[i].contains(patient)) {
                List<Integer> route = genes[i];
                int patientPosition = route.indexOf(patient);

                int j = routeEndPoint[i];
                while(routeEndPoint[i] != patientPosition && j < route.size()) {
                    int curPatient;
                    if(j==0)
                        curPatient = -1;
                    else
                        curPatient = route.get(j-1);
                    if(patientIsAssigned(genes,i,curPatient,route.get(j),totalTravelCost,routesCurrentTime,highestAndTotalTardiness,track)==Double.POSITIVE_INFINITY)
                        return Integer.MAX_VALUE;
                    j++;
                }
                return i;
            }
        }
        return Integer.MAX_VALUE;
    }

    private boolean SwapRoutes(Patient p, int route1, int route2) {
        Required_Caregiver[] requiredCaregivers = p.getRequired_caregivers();
        String service1 = requiredCaregivers[0].getService();
        String service2 = requiredCaregivers[1].getService();

        Set<Integer> service1Routes = dataset.getQualifiedCaregiver(service1);
        Set<Integer> service2Routes = dataset.getQualifiedCaregiver(service2);

        return (service1Routes.contains(route2) && !service2Routes.contains(route2)) ||
                (service2Routes.contains(route1) && !service1Routes.contains(route1)) ||
                (service1Routes.contains(route2) && !service1Routes.contains(route1));
    }

    private boolean noEvaluationConflicts(List<Integer> c1Route, List<Integer> c2Route, int m, int n) {
        return conflictCheck(c1Route, c2Route, m, n);
    }
}
