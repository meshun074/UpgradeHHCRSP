package org.example.GA;

import org.example.Data.InstancesClass;
import org.example.Data.Patient;
import org.example.Data.Required_Caregiver;

import java.util.*;

import static org.example.GA.EvaluationFunction.removeAffectedPatient;
import static org.example.GA.GeneticAlgorithm.conflictCheck;

public class LocalSearchNew implements Runnable {
    @Override
    public void run() {
        ga.getLSChromosomes().put(r, search());
    }

    private GeneticAlgorithm ga;
    private Chromosome ch;
    private final Random rand;
    private final int r;
    private final int gen;
    private static Patient[] allPatients;
    private static InstancesClass dataset;
    private static int allCaregivers;
    private static double[][] distances;
    private final int[] routeEndPoint;
    private final double[] routesCurrentTime;
    private final double[] highestAndTotalTardiness;
    private final Set<Integer> track;

    public LocalSearchNew(GeneticAlgorithm ga, Chromosome ch, int r, Random rand, int gen) {
        this.ga = ga;
        this.ch = ch;
        this.rand = rand;
        this.r = r;
        this.gen = gen;
        this.routeEndPoint = new int[allCaregivers];
        this.routesCurrentTime = new double[allCaregivers];
        this.highestAndTotalTardiness = new double[2];
        this.track = new HashSet<>(100);
    }

    public static void initialize(InstancesClass data) {
        allPatients = data.getPatients();
        distances = data.getDistances();
        allCaregivers = data.getCaregivers().length;
        dataset = data;
    }

    public Chromosome search() {
        Chromosome best = ch;
        boolean continueLS;
        boolean genCont = gen >= 100;
        int counter = 0;
        do {
            continueLS = false;
            ch = localSearch(ch);
            if (best.getFitness() - ch.getFitness() > 0.001 || ch.getFitness() <= best.getFitness() && rand.nextBoolean()) {
                best = ch;
                continueLS = true;
            } else {
                if (genCont && counter < 6) {
                    continueLS = true;
                }
            }
            counter++;
        } while (continueLS);
        return best;
    }

    private Chromosome localSearch(Chromosome ch) {
        List<Integer>[] p1Routes, c1Routes;
        int patientLength = allPatients.length;
        int size = (int) (patientLength * (0.4+0.1*Math.random()));
        Set<Integer> selectRoute = new LinkedHashSet<>(size);
        Random rand = new Random();
        while (selectRoute.size() < size) {
            int sp = rand.nextInt(patientLength);
            selectRoute.add(sp);
        }
        p1Routes = ch.getGenes();
        c1Routes = new ArrayList[p1Routes.length];
        //Removing patients of selected route from parent routes
        for (int i = 0; i < p1Routes.length; i++) {
            List<Integer> route1 = new ArrayList<>(p1Routes[i].size());
            for (int j = 0; j < p1Routes[i].size(); j++) {
                int patient = p1Routes[i].get(j);
                if (!selectRoute.contains(patient)) {
                    route1.add(patient);
                }
            }
            c1Routes[i] = route1;
        }
        //Inserting removed route
        List<Integer> route2 = new ArrayList<>(selectRoute);
        Chromosome cTemp = new Chromosome(c1Routes, 0.0, true);
        EvaluationFunction.Evaluate(cTemp);
        boolean isInvalid;
        for (int i = 0; i < route2.size(); i++) {
            isInvalid = cTemp.getFitness() == Double.POSITIVE_INFINITY;
            int patient = route2.get(i);
            Patient p = allPatients[patient];
            double bestCost = Double.MAX_VALUE;
            int bestFirst = -1;
            int bestSecond = -1;
            int bestM = -1;
            int bestN = -1;
            cTemp.buildPatientRouteMap();
            Shift[] shifts = cTemp.getCaregiversRouteUp();
            if (p.getRequired_caregivers().length > 1) {
                List<CaregiverPair> caregiverPairs = p.getAllPossibleCaregiverCombinations();
                for (int x = 0; x < caregiverPairs.size(); x++) {
                    CaregiverPair caregiverPair = caregiverPairs.get(x);
                    int first = caregiverPair.getFirst();
                    int second = caregiverPair.getSecond();
                    for (int m = 0; m <= c1Routes[first].size(); m++) {
                        for (int n = 0; n <= c1Routes[second].size(); n++) {
                            if (noEvaluationConflicts(c1Routes[first], c1Routes[second], m, n)) {
                                double tempCost = calMoveCost(first, m, second, n, patient, cTemp, bestCost, shifts, isInvalid);
                                if (bestCost == Double.MAX_VALUE || bestCost - tempCost > 0.001 || tempCost <= bestCost && rand.nextBoolean()) {
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
                    isInvalid = cTemp.getFitness() == Double.POSITIVE_INFINITY;
                    cTemp = swap(cTemp, isInvalid, patient, bestFirst, bestSecond, bestM, bestN);
                }
            } else {
                List<CaregiverPair> caregiverPairs = p.getAllPossibleCaregiverCombinations();
                for (int x = 0; x < caregiverPairs.size(); x++) {
                    CaregiverPair caregiverPair = caregiverPairs.get(x);
                    int first = caregiverPair.getFirst();
                    for (int k = 0; k <= c1Routes[first].size(); k++) {
                        double tempCost = calMoveCost(first, k, -1, -1, patient, cTemp, bestCost, shifts, isInvalid);
                        if (bestCost == Double.MAX_VALUE || bestCost - tempCost > 0.001 || tempCost <= bestCost && rand.nextBoolean()) {
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
                    isInvalid = cTemp.getFitness() == Double.POSITIVE_INFINITY;
                    cTemp = swap(cTemp, isInvalid, patient, bestFirst, -1, bestM, -1);
                }
            }
        }
        return cTemp;
    }

    private double calMoveCost(int first, int m, int second, int n, int patient, Chromosome cTemp, double bestCost, Shift[] shifts, boolean isInvalid) {
        Arrays.fill(routeEndPoint, -1);
        if (isInvalid) {
            List<Integer>[] c1Routes = cTemp.getGenes();
            c1Routes[first].add(m, patient);
            if (second != -1) {
                c1Routes[second].add(n, patient);
            }
            Chromosome temp = new Chromosome(c1Routes, 0.0, true);
            c1Routes[first].remove(Integer.valueOf(patient));
            if (second != -1) {
                c1Routes[second].remove(Integer.valueOf(patient));
            }
            EvaluationFunction.Evaluate(temp);
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
        if (second != -1) {
            genes[second].add(n, patient);
        }
        for (int i = 0; i < routeEndPoint.length; i++) {
            List<Integer> route = genes[i];
            int routeEnd = routeEndPoint[i];
            if (i == first || i == second) {
                for (int j = routeEnd; j <= route.size(); j++) {
                    if (j == 0) {
                        int nextIndex = route.get(j) + 1;
                        totalTravelCost += distances[0][nextIndex];
                    } else if (j == route.size()) {
                        int prevIndex = route.get(j - 1) + 1;
                        totalTravelCost += distances[prevIndex][0];
                    } else {
                        int nextIndex = route.get(j) + 1;
                        int prevIndex = route.get(j - 1) + 1;
                        totalTravelCost += distances[prevIndex][nextIndex];
                    }

                }
            }
        }
        totalTravelCost = 1 / 3d * totalTravelCost;

        double solutionCost = totalTravelCost + (1 / 3d * highestAndTotalTardiness[0]) + (1 / 3d * highestAndTotalTardiness[1]);
        if (solutionCost > bestCost) {
            genes[first].remove(Integer.valueOf(patient));
            if (second != -1) {
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
                        solutionCost = patientIsAssigned(genes, i, route.get(j - 1), route.get(j), totalTravelCost, routesCurrentTime, highestAndTotalTardiness, track);
                    }
                    if (solutionCost == Double.POSITIVE_INFINITY || solutionCost > bestCost) {
                        genes[first].remove(Integer.valueOf(patient));
                        if (second != -1) {
                            genes[second].remove(Integer.valueOf(patient));
                        }
                        return solutionCost;
                    }

                    track.clear();
                }
            }
        }
        genes[first].remove(Integer.valueOf(patient));
        if (second != -1) {
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

            int route2 = findOtherCaregiver(patient, route1, genes, totalTravelCost, routesCurrentTime, highestAndTotalTardiness, track);
            if (route2 > allCaregivers - 1) {
                return Double.POSITIVE_INFINITY;
            }

            Required_Caregiver[] requiredCaregivers = p.getRequired_caregivers();

            int position2 = routeEndPoint[route2] - 1;
            int curPatient2;
            if (position2 == -1) {
                curPatient2 = -1;
            } else {
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
            } else {
                double startTime = Math.max(startTime1, startTime2);
                double tardiness = Math.max(0, startTime - timeWindow[1]);
                highestAndTotalTardiness[0] = Math.max(highestAndTotalTardiness[0], tardiness);
                highestAndTotalTardiness[1] += (2 * tardiness);
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
                while (routeEndPoint[i] != patientPosition && j < route.size()) {
                    int curPatient;
                    if (j == 0)
                        curPatient = -1;
                    else
                        curPatient = route.get(j - 1);
                    if (patientIsAssigned(genes, i, curPatient, route.get(j), totalTravelCost, routesCurrentTime, highestAndTotalTardiness, track) == Double.POSITIVE_INFINITY)
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

    private Chromosome swap(Chromosome ch, boolean isInvalid, int patient, int first, int second, int firstPosition, int secondPosition) {
        ch.buildPatientRouteMap();
        double bestCost = Double.MAX_VALUE;
        double currentCost = ch.getFitness();
        List<Integer>[] routes = ch.getGenes();
        Shift[] shifts = ch.getCaregiversRouteUp();

        if (second != -1) {
            List<Integer> route1 = routes[first];
            List<Integer> route2 = routes[second];
            int bestZ = -1;
            int bestFirstPatient = -1;
            int bestSecondPatient = -1;
            int bestL = -1;
            for (int z = 0; z < route1.size(); z++) {
                if (Math.abs(z - firstPosition) > 1) {
                    int p3 = route1.get(z);
                    double tempCost = calSwapMoveCost(first, firstPosition, -1, -1, patient, p3, z, -1, -1, ch, bestCost, shifts, isInvalid);
                    if (bestCost == Double.MAX_VALUE || bestCost - tempCost > 0.001 || tempCost <= bestCost && rand.nextBoolean()) {
                        bestCost = tempCost;
                        bestZ = z;
                        bestFirstPatient = p3;
                    }
                }
            }

            for (int l = 0; l < route2.size(); l++) {
                if (Math.abs(l - secondPosition) > 1) {
                    int p4 = route2.get(l);
                    double tempCost = calSwapMoveCost(-1, -1, second, secondPosition, patient, -1, -1, p4, l, ch, bestCost, shifts, isInvalid);
                    if (bestCost == Double.MAX_VALUE || bestCost - tempCost > 0.001 || tempCost <= bestCost && rand.nextBoolean()) {
                        bestCost = tempCost;
                        bestL = l;
                        bestSecondPatient = p4;
                        bestZ = -1;
                        bestFirstPatient = -1;
                    }
                }
            }
            for (int z = 0; z < route1.size(); z++) {
                if (Math.abs(z - firstPosition) > 1) {
                    for (int l = 0; l < route2.size(); l++) {
                        if (Math.abs(l - secondPosition) > 1) {
                            int p3 = route1.get(z);
                            int p4 = route2.get(l);
                            double tempCost = calSwapMoveCost(first, firstPosition, second, secondPosition, patient, p3, z, p4, l, ch, bestCost, shifts, isInvalid);
                            if (bestCost == Double.MAX_VALUE || bestCost - tempCost > 0.001 || tempCost <= bestCost && rand.nextBoolean()) {
                                bestCost = tempCost;
                                bestL = l;
                                bestSecondPatient = p4;
                                bestZ = z;
                                bestFirstPatient = p3;
                            }
                        }
                    }
                }
            }
            if (bestCost != Double.MAX_VALUE && currentCost - bestCost > 0.001 || bestCost <= currentCost && bestCost != Double.MAX_VALUE && rand.nextBoolean()) {
                if (bestZ != -1) {
                    routes[first].set(firstPosition, bestFirstPatient);
                    routes[first].set(bestZ, patient);
                }
                if (bestL != -1) {
                    routes[second].set(secondPosition, bestSecondPatient);
                    routes[second].set(bestL, patient);
                }
                EvaluationFunction.Evaluate(ch);
            }
        } else {
            List<Integer> route1 = routes[first];
            int bestI = -1;
            int bestPatient = -1;
            for (int i = 0; i < route1.size(); i++) {
                if (Math.abs(i - firstPosition) > 1) {
                    int p2 = route1.get(i);

                    double tempCost = calSwapMoveCost(first, firstPosition, -1, -1, patient, p2, i, -1, -1, ch, bestCost, shifts, isInvalid);
                    if (bestCost == Double.MAX_VALUE || bestCost - tempCost > 0.001 || tempCost <= bestCost && rand.nextBoolean()) {
                        bestCost = tempCost;
                        bestI = i;
                        bestPatient = p2;
                    }
                }
            }
            if (bestCost != Double.MAX_VALUE && currentCost - bestCost > 0.001 || bestCost <= currentCost && bestCost != Double.MAX_VALUE && rand.nextBoolean()) {
                routes[first].set(firstPosition, bestPatient);
                routes[first].set(bestI, patient);
                EvaluationFunction.Evaluate(ch);
            }
        }
        return ch;
    }

    private double calSwapMoveCost(int first, int m, int second, int n, int patient, int patient1, int patient1Position, int patient2, int patient2Position, Chromosome cTemp, double bestCost, Shift[] shifts, boolean isInvalid) {
        Arrays.fill(routeEndPoint, -1);
        if (isInvalid) {
            List<Integer>[] c1Routes = cTemp.getGenes();
            if (patient1 != -1) {
                c1Routes[first].set(m, patient1);
                c1Routes[first].set(patient1Position, patient);
            }
            if (patient2 != -1) {
                c1Routes[second].set(n, patient2);
                c1Routes[second].set(patient2Position, patient);
            }
            Chromosome temp = new Chromosome(c1Routes, 0.0, false);
            EvaluationFunction.Evaluate(temp);
            System.out.println(temp + " Solution cost is " + temp.getFitness() + " travel cost is " + temp.getFitness() + " tardiness " + temp.getTotalTardiness() + " highest tardiness " + temp.getHighestTardiness());
            if (patient1 != -1) {
                c1Routes[first].set(m, patient);
                c1Routes[first].set(patient1Position, patient1);
            }
            if (patient2 != -1) {
                c1Routes[second].set(n, patient);
                c1Routes[second].set(patient2Position, patient2);
            }
            return temp.getFitness();
        }
        int size = 0;
        int newFirstPosition = Math.min(m, patient1Position);
        if (patient1 != -1) {
            routeEndPoint[first] = newFirstPosition;
            size++;
        }
        int newSecondPosition = Math.min(n, patient2Position);
        if (patient2 != -1) {
            routeEndPoint[second] = newSecondPosition;
            size++;
        }

        int[] routeMove = new int[size];
        int[] positionMove = new int[size];
        int counter = 0;
        if (patient1 != -1) {
            routeMove[counter] = first;
            positionMove[counter] = newFirstPosition;
            counter++;
        }
        if (patient2 != -1) {
            routeMove[counter] = second;
            positionMove[counter] = newSecondPosition;
        }
//        System.out.println("patients:"+patient1+","+patient2+","+counter);
//        for (int i = 0; i < routeMove.length; i++) {
//            System.out.println("Route: "+routeMove[i]);
//        }
//        for (int i = 0; i < positionMove.length; i++) {
//            System.out.println("position: "+positionMove[i]);
//        }
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
            if (i == first && patient1 != -1 || i == second && patient2 != -1) {
                int index = routeEndPoint[i];
                totalTravelCost += travelCost.get(index);
//                System.out.println("sample travel: "+travelCost.get(index));
            } else {
                totalTravelCost += travelCost.get(travelCost.size() - 1);
            }
        }

//        System.out.println("After total travel cost: "+totalTravelCost);

        //Distance calculation
        List<Integer>[] genes = cTemp.getGenes();
        if (patient1 != -1) {
            genes[first].set(m, patient1);
            genes[first].set(patient1Position, patient);
        }
        if (patient2 != -1) {
            genes[second].set(n, patient2);
            genes[second].set(patient2Position, patient);
        }
        for (int i = 0; i < routeEndPoint.length; i++) {
            List<Integer> route = genes[i];
            int routeEnd = routeEndPoint[i];
            if (i == first && patient1 != -1 || i == second && patient2 != -1) {
                for (int j = routeEnd; j <= route.size(); j++) {
                    if (j == 0) {
                        int nextIndex = route.get(j) + 1;
                        totalTravelCost += distances[0][nextIndex];
                    } else if (j == route.size()) {
                        int prevIndex = route.get(j - 1) + 1;
                        totalTravelCost += distances[prevIndex][0];
                    } else {
                        int nextIndex = route.get(j) + 1;
                        int prevIndex = route.get(j - 1) + 1;
                        totalTravelCost += distances[prevIndex][nextIndex];
                    }

                }
            }
        }
        totalTravelCost = 1 / 3d * totalTravelCost;

        double solutionCost = totalTravelCost + (1 / 3d * highestAndTotalTardiness[0]) + (1 / 3d * highestAndTotalTardiness[1]);
        if (solutionCost > bestCost) {
//            System.out.println(cTemp+" Solution cost is " + solutionCost+ " travel cost is " + totalTravelCost*3+" tardiness "+ highestAndTotalTardiness[1]+ " highest tardiness "+ highestAndTotalTardiness[0]);
            if (patient1 != -1) {
                genes[first].set(m, patient);
                genes[first].set(patient1Position, patient1);
            }
            if (patient2 != -1) {
                genes[second].set(n, patient);
                genes[second].set(patient2Position, patient2);
            }

            return solutionCost;
        }
//        Set<Integer> track = new HashSet<>(100);

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
                        solutionCost = patientIsAssigned(genes, i, route.get(j - 1), route.get(j), totalTravelCost, routesCurrentTime, highestAndTotalTardiness, track);
                    }
                    if (solutionCost == Double.POSITIVE_INFINITY || solutionCost > bestCost) {
//                        System.out.println(cTemp+" Solution cost is " + solutionCost+ " travel cost is " + totalTravelCost*3+" tardiness "+ highestAndTotalTardiness[1]+ " highest tardiness "+ highestAndTotalTardiness[0]);
                        if (patient1 != -1) {
                            genes[first].set(m, patient);
                            genes[first].set(patient1Position, patient1);
                        }
                        if (patient2 != -1) {
                            genes[second].set(n, patient);
                            genes[second].set(patient2Position, patient2);
                        }

                        return solutionCost;
                    }

                    track.clear();
                }
            }
        }
//        System.out.println(cTemp+" Solution cost is " + solutionCost+ " travel cost is " + totalTravelCost*3+" tardiness "+ highestAndTotalTardiness[1]+ " highest tardiness "+ highestAndTotalTardiness[0]);
        if (patient1 != -1) {
            genes[first].set(m, patient);
            genes[first].set(patient1Position, patient1);
        }
        if (patient2 != -1) {
            genes[second].set(n, patient);
            genes[second].set(patient2Position, patient2);
        }
//        System.out.println("travel cost: " + totalTravelCost*3 +" highest Tardiness: " + highestAndTotalTardiness[0] + " tardiness: " + highestAndTotalTardiness[1]);

        return solutionCost;
    }

    private boolean noEvaluationConflicts(List<Integer> c1Route, List<Integer> c2Route, int m, int n) {
        return conflictCheck(c1Route, c2Route, m, n);
    }
}
