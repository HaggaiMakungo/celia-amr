const state = {
  candidates: [],
  selectedCandidate: null,
  selectedScenario: "Targeted response",
  assessment: null,
};

const ANALYSIS_ANIMATION_MS = 760;
const SCRAMBLE_CHARACTERS = "ABCDEFGHIJKLMNOPQRSTUVWXYZ0123456789";
const liveBackend = new URLSearchParams(window.location.search).get("live") === "1";
const demoData = window.CELIA_DEMO_DATA;

const pages = [...document.querySelectorAll(".page")];
const navItems = [...document.querySelectorAll(".nav-item")];

function openPage(view) {
  pages.forEach((page) => page.classList.toggle("is-active", page.id === view));
  navItems.forEach((item) => {
    const active = item.dataset.view === view;
    item.classList.toggle("is-active", active);
    if (active) item.setAttribute("aria-current", "page");
    else item.removeAttribute("aria-current");
  });
  history.replaceState(null, "", `#${view}`);
  document.querySelector("#main-content").focus?.();
}

document.querySelectorAll("[data-view]").forEach((control) => {
  control.addEventListener("click", () => openPage(control.dataset.view));
});

function renderCandidates() {
  const target = document.querySelector("#candidate-list");
  target.innerHTML = state.candidates.map((candidate) => `
    <button class="candidate-option ${candidate === state.selectedCandidate ? "is-selected" : ""}" type="button" data-candidate="${candidate.name}">
      <span><b>${candidate.name}</b><small>${candidate.status} · ${candidate.fit} fit</small></span><em>${candidate.score}</em>
    </button>`).join("");

  target.querySelectorAll("[data-candidate]").forEach((button) => {
    button.addEventListener("click", () => {
      state.selectedCandidate = state.candidates.find((candidate) => candidate.name === button.dataset.candidate);
      updateCandidateCard();
      renderCandidates();
    });
  });
}

function updateCandidateCard() {
  const candidate = state.selectedCandidate;
  if (!candidate) return;
  document.querySelector("#candidate-name").textContent = candidate.name;
  document.querySelector("#candidate-status").textContent = candidate.status.toUpperCase();
  document.querySelector("#candidate-score").textContent = candidate.score;
  document.querySelector("#candidate-fit").textContent = candidate.fit;
  document.querySelector("#candidate-copy").textContent = candidate.name === "ZMB-041"
    ? "High-confidence molecular fit with a tractable next validation."
    : "A viable alternate candidate for additional evidence review.";
  document.querySelector("#brief-candidate").textContent = `${candidate.name} selected for follow-up`;
}

function resetDemo() {
  state.assessment = null;
  state.selectedScenario = "Targeted response";
  document.querySelector("#analysis-form").reset();

  const result = document.querySelector("#analysis-result");
  result.hidden = true;
  result.innerHTML = "";
  document.querySelector("#brief-genomic").textContent = "Run an assessment to add genomic evidence.";
  document.querySelector("#recommendation-copy").textContent = "Carry the targeted response into the research brief.";
  document.querySelector("#brief-risk").textContent = "Targeted response preferred";
  document.querySelectorAll(".scenario").forEach((button) => {
    button.classList.toggle("is-selected", button.dataset.scenario === "Targeted response");
  });

  if (state.candidates.length) {
    state.selectedCandidate = state.candidates[0];
    updateCandidateCard();
    renderCandidates();
  }
  openPage("projects");
}

function wait(milliseconds) {
  return new Promise((resolve) => window.setTimeout(resolve, milliseconds));
}

function createAssessmentPreview(gene) {
  const result = document.querySelector("#analysis-result");
  result.hidden = false;
  result.innerHTML = `
    <p class="eyebrow">ASSESSMENT · ${gene.toUpperCase()}</p>
    <div class="result-summary result-summary--preview"><span class="result-score result-score--rolling" aria-hidden="true">00</span><div><b class="preview-label">CALCULATING</b><p class="result-meta preview-copy">Preparing assessment</p></div></div>
    <div class="analysis-preview" aria-live="polite">
      <span class="analysis-preview__label">BUILDING RESULT</span>
      <strong class="analysis-preview__text">Preparing assessment</strong>
      <i class="analysis-preview__cursor" aria-hidden="true"></i>
      <div class="analysis-preview__lines" aria-hidden="true"><i></i><i></i><i></i></div>
    </div>`;

  const score = result.querySelector(".result-score--rolling");
  const text = result.querySelector(".analysis-preview__text");
  const messages = ["Mapping evidence", "Resolving signal", "Preparing assessment"];
  let tick = 0;
  const timer = window.setInterval(() => {
    score.textContent = String(Math.floor(Math.random() * 100)).padStart(2, "0");
    const message = messages[Math.floor(tick / 4) % messages.length];
    const revealed = tick % (message.length + 1);
    text.textContent = [...message].map((character, index) => {
      if (index < revealed || character === " ") return character;
      return SCRAMBLE_CHARACTERS[Math.floor(Math.random() * SCRAMBLE_CHARACTERS.length)];
    }).join("");
    tick += 1;
  }, 75);

  return () => window.clearInterval(timer);
}

function animateScore(scoreElement, target) {
  const startTime = performance.now();
  const duration = 620;
  const update = (now) => {
    const progress = Math.min((now - startTime) / duration, 1);
    const eased = 1 - Math.pow(1 - progress, 3);
    scoreElement.textContent = Math.round(target * eased);
    if (progress < 1) window.requestAnimationFrame(update);
  };
  window.requestAnimationFrame(update);
}

function renderAssessment(payload) {
  const { assessment, public_health_context: context } = payload;
  state.assessment = assessment;
  const mutations = assessment.mutation_report.filter((row) => row.is_driving_mutation || row.status === "resistant");
  const result = document.querySelector("#analysis-result");
  result.hidden = false;
  result.innerHTML = `
    <p class="eyebrow">ASSESSMENT · ${assessment.gene.toUpperCase()}</p>
    <div class="result-summary analysis-reveal" style="--reveal-delay: 0ms"><span class="result-score" id="result-score">0</span><div><b>${assessment.band}</b><p class="result-meta">${assessment.drug} · genomic evidence ${assessment.genomic_component}/100</p></div></div>
    <h2 class="analysis-reveal" style="--reveal-delay: 90ms">${assessment.notes[0] || "Genomic evidence ready"}</h2>
    <ul class="mutation-list analysis-reveal" style="--reveal-delay: 180ms">${mutations.map((row) => `<li><span><b>${row.mutation}</b><br><small>${row.status.replaceAll("_", " ")}</small></span><small>WHO grade ${row.who_confidence_grade}</small></li>`).join("")}</ul>
    <p class="context-note analysis-reveal" style="--reveal-delay: 270ms">Zambia context, ${context.year}: ${Number(context.incident_cases).toLocaleString()} estimated incident TB cases. This context is shown separately and does not alter the sample resistance score.</p>`;
  animateScore(result.querySelector("#result-score"), assessment.resistance_score);
  document.querySelector("#brief-genomic").textContent = `${assessment.gene} ${assessment.notes[0].replace("Driving mutation ", "")} detected`;
}

async function runAssessment(event) {
  event.preventDefault();
  const button = document.querySelector("#run-analysis");
  const result = document.querySelector("#analysis-result");
  const gene = new FormData(event.currentTarget).get("gene");
  button.disabled = true;
  button.innerHTML = "Running assessment…";
  const stopPreview = createAssessmentPreview(gene);
  const startedAt = performance.now();
  try {
    let payload;
    if (liveBackend) {
      const response = await fetch("/api/analysis", {
        method: "POST",
        headers: { "Content-Type": "application/json" },
        body: JSON.stringify({ gene }),
      });
      payload = await response.json();
      if (!response.ok) throw new Error(payload.error || "The assessment could not be completed.");
    } else {
      payload = demoData.assessments[gene];
      if (!payload) throw new Error("The selected demo result is unavailable.");
    }
    await wait(Math.max(0, ANALYSIS_ANIMATION_MS - (performance.now() - startedAt)));
    stopPreview();
    renderAssessment(payload);
  } catch (error) {
    stopPreview();
    result.innerHTML = `<p class="eyebrow">ASSESSMENT</p><h2>Assessment unavailable.</h2><p>${error.message}</p>`;
  } finally {
    button.disabled = false;
    button.innerHTML = "Run assessment <b>→</b>";
  }
}

document.querySelector("#analysis-form").addEventListener("submit", runAssessment);
document.querySelector("#reset-demo").addEventListener("click", resetDemo);

document.querySelectorAll(".scenario").forEach((button) => {
  button.addEventListener("click", () => {
    document.querySelectorAll(".scenario").forEach((item) => item.classList.toggle("is-selected", item === button));
    state.selectedScenario = button.dataset.scenario;
    document.querySelector("#recommendation-copy").textContent = state.selectedScenario === "Targeted response"
      ? "Carry the targeted response into the research brief."
      : "Review the higher-risk baseline with the programme team.";
    document.querySelector("#brief-risk").textContent = `${state.selectedScenario} preferred`;
  });
});

document.querySelector("#print-brief").addEventListener("click", () => window.print());

async function initialise() {
  try {
    let candidates;
    if (liveBackend) {
      const response = await fetch("/api/project");
      const payload = await response.json();
      if (!response.ok) throw new Error("Project data is unavailable.");
      candidates = payload.candidates;
    } else {
      candidates = demoData.candidates;
    }
    state.candidates = candidates;
    state.selectedCandidate = state.candidates[0];
    updateCandidateCard();
    renderCandidates();
  } catch {
    document.querySelector("#candidate-list").innerHTML = "<p>Project data is unavailable.</p>";
  }
  const view = window.location.hash.slice(1);
  if (pages.some((page) => page.id === view)) openPage(view);
}

initialise();
