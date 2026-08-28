/* Static values for the public GitHub Pages demonstration. */
window.CELIA_DEMO_DATA = {
  candidates: [
    { name: "ZMB-041", score: "−9.2", fit: "High", status: "Lead candidate" },
    { name: "ZMB-117", score: "−8.7", fit: "Promising", status: "Evidence review" },
    { name: "ZMB-203", score: "−8.4", fit: "Promising", status: "Evidence review" },
  ],
  assessments: {
    katG: {
      assessment: {
        gene: "katG",
        resistance_score: 95,
        band: "High confidence resistant",
        drug: "Isoniazid",
        genomic_component: 95,
        notes: ["Sample assessment complete"],
        mutation_report: [
          { mutation: "S315T", status: "resistant", is_driving_mutation: true, who_confidence_grade: 1 },
          { mutation: "D419A", status: "resistant", is_driving_mutation: false, who_confidence_grade: 4 },
        ],
      },
      public_health_context: { year: 2023, incident_cases: 59000 },
    },
    rpoB: {
      assessment: {
        gene: "rpoB",
        resistance_score: 92,
        band: "High confidence resistant",
        drug: "Rifampicin",
        genomic_component: 92,
        notes: ["Sample assessment complete"],
        mutation_report: [
          { mutation: "S450L", status: "resistant", is_driving_mutation: true, who_confidence_grade: 1 },
          { mutation: "D435V", status: "resistant", is_driving_mutation: false, who_confidence_grade: 2 },
        ],
      },
      public_health_context: { year: 2023, incident_cases: 59000 },
    },
  },
};
