# Спектральная щель нормализованного лапласиана в разреженных случайных графах с тяжёлохвостыми степенями: асимптотика, bottleneck-структуры и связь с модульностью биологических сетей

> ## Статус: проект закрыт — содержательной новизны нет
>
> В ходе работы эмпирически воспроизведён и численно детализирован **известный**
> результат: **Samukhin, Dorogovtsev, Mendes**, *Laplacian spectra of complex
> networks and random walks on them*, **PRE 77, 036115 (2008)** — нижний край
> спектра лапласиана разреженного случайного графа определяется **минимальной
> степенью** (`λ_c > 0` при `q_m > 2`, `λ_c = 0` при `q_m ≤ 2`), а тяжёлый хвост
> вторичен. Конечный размер тоже известен: `λ₂(N) ∼ (ln N)^{−2}` при `q_m ≤ 2`
> и `λ₂(N) − λ_c ∼ (ln ln N)^{−2}` при `q_m > 2`.
>
> **Нового результата в этой нише получить нельзя.** Репозиторий сохранён как
> воспроизводимый **инструмент/бенчмарк**: генераторы, спектральная диагностика,
> ER/toy-контроли, тесты. Детали и разбор литературы — в `LITERATURE.md` и
> `CHANGELOG.md`.

Рабочий документ проекта. Постановка, формализация, диагностики, алгоритмы, план и открытые вопросы.

Обозначения: λ₂, λ₃ — собственные значения нормализованного лапласиана; φ — conductance (проводимость); ρ — доля объёма; Δ — параметр конкуренции bottleneck'ов; d̄ — средняя степень; γ — показатель степенного хвоста; d_min, d_max — минимальная и максимальная степени; G_giant — гигантская компонента; G_2core — 2-ядро (2-core).

---

## 1. Мотивация и исходная постановка

**Исходная идея:** получить точную асимптотику второй собственной величины λ₂ нормализованного лапласиана для разреженных случайных графов с тяжёлохвостым распределением степеней P(D = k) ∼ C·k^(−γ) в режиме d̄_min = O(1), не покрытом теоремой Chung–Lu–Vu.

**Проблема исходной постановки:** она слишком широкая и в части параметров тривиально неверна. В разреженном режиме λ₂ может быть нулём из-за несвязности, а после удаления мелких компонент — иметь совершенно другой порядок. Это не техническая мелочь, а центральная часть исследования.

**Корректировка:** вместо одной формулы λ₂ ∼ f(γ, d̄, n) нужно строить **фазовую диаграмму порядка λ₂** в зависимости от:

    γ,   d_min,   d_max,   d̄,   n,

и, что особенно важно, **структуры низкостепенного хвоста**.

---

## 2. Что уже известно из литературы

- **Chung–Lu–Vu (2003):** для случайных графов с заданными ожидаемыми степенями при d̄_min ≫ ln² n спектр нормализованного лапласиана концентрируется около предсказываемого поведения; результаты включают power-law графы. Оценка щели:

      λ(G) ≥ 1 − 4·d̄^(−1/2) − d̄_min^(−1)·ln² n.

  При d̄_min ≤ ln² n эта оценка становится бессмысленной.

- **Coja-Oghlan–Lanka (2009):** для разреженного режима d̄_min = O(1) показано, что существует большой core, у которого spectral gap нормализованного лапласиана остаётся отделённой от нуля:

      λ₂(L_core) ≥ 1 − c₀·d̄_min^(−1/2)

  с высокой вероятностью. Речь идёт именно о **core**, а не обо всём исходном графе.

- **Для G(n, p) при ограниченной средней степени d:** spectral gap всего графа нормализованного лапласиана имеет λ₂ = o(1), хотя giant/core может иметь gap порядка единицы.

**Ключевой вывод:**

> **λ₂(G) ≁ λ₂(core(G)).**

---

## 3. Почему одного γ недостаточно

Пусть P(D = k) ∼ C·k^(−γ). Это определяет хвост при больших k. Но λ₂ очень чувствительна к **маленьким степеням**.

**Контрпример:** две последовательности с одинаковым γ = 2.5:

- Модель A: P(D = 1) > 0.
- Модель B: P(D = 1) = 0, P(D ≥ 2) = 1.

У них одинаковый асимптотический хвост, но совершенно разное поведение λ₂. В первой модели появляются листья и конечные компоненты; во второй граф может оказаться существенно более связным.

> **γ не определяет λ₂.**

**Жёсткий контрпример:** конфигурационная модель с P(D = 1) = p₁ > 0. При ограниченной средней степени такие графы естественно содержат мелкие компоненты и dangling structures. Если граф несвязен, то λ₂(L) = 0. Следовательно, никакая нетривиальная формула вида λ₂ ∼ f(γ) > 0 для всего графа невозможна.

---

## 4. Три спектральные щели

| Объект | Обозначение | Что характеризует |
|---|---|---|
| Полный граф | λ₂(G) | Глобальная связность; λ₂ = 0 при несвязности |
| Giant component | λ₂(G_giant) | Содержательная величина для sparse heavy-tailed сетей |
| 2-core | λ₂(G_2core) | Ненулевая asymptotic lower bound даже при d̄_min = O(1) |

**Гипотеза о порядке:**

> **λ₂(G) ≪ λ₂(G_giant) ≲ λ₂(G_2core).**

---

## 5. Три механизма возникновения малой λ₂

1. **Мелкие компоненты:** λ₂(G) = 0.
2. **Dangling trees:** длинные «усы» создают малый conductance; по Чигеру λ₂ ≤ 2φ, поэтому достаточно длинный тонкий отросток может сделать λ₂ → 0.
3. **Настоящая модульность:** две большие dense communities с малым cut дают φ ≪ 1 и λ₂ ≪ 1.

**Ключевое разделение:**

    малое λ₂  ⟹  существует слабый глобальный cut

но

    малое λ₂  ⇏  биологическая модульность.

Нужно отделить **community bottleneck** от **peripheral bottleneck**.

---

## 6. Формализация двух типов bottleneck

Пусть

    vol(S) = Σ_{v ∈ S} d_v,        φ(S) = |∂S| / vol(S).

**Периферийный bottleneck:**

    φ_periph(G) = min_{ S ⊂ V : 0 < vol(S) ≤ ε·vol(V) } φ(S),

где ε ≪ 1, например 0.05.

**Модульный bottleneck:**

    φ_mod(G) = min_{ S ⊂ V : ε·vol(V) ≤ vol(S) ≤ ½·vol(V) } φ(S).

Тогда

    φ(G) = min{ φ_periph, φ_mod }.

Через Чигера:

    φ(G)² / 2  ≤  λ₂  ≤  2·φ(G).

То есть мы получаем не «разложение λ₂», а **классификацию механизма, создающего глобальный bottleneck**.

---

## 7. Диагностические параметры

**Параметр конкуренции:**

> **Δ(G) = φ_mod(G) / φ_periph(G)**

с интерпретацией:

- Δ ≫ 1 — периферийный bottleneck дешевле;
- Δ ≪ 1 — модульный bottleneck дешевле;
- Δ ≈ 1 — переходный режим.

**Для устойчивости численных оценок использовать:**

    log Δ = log φ_mod − log φ_periph.

**Размер минимального спектрального cut:**

    ρ(S) = vol(S) / vol(V),        ρ_* = ρ(S_*),

где S_* = argmin_S φ(S).

Интерпретация:

- ρ_* = O(n^(−1)) — маленькая локальная структура;
- ρ_* = O(1) — макроскопическое разделение графа.

**Диагностическая плоскость:** (log Δ, ρ_*).

**Спектральный cut vs комбинаторный:**

    ρ₂ = vol({ v : f₂(v) > 0 }) / vol(V).

Если ρ₂ ≈ ρ_* — спектральный cut совпадает с комбинаторным. Если нет — λ₂ контролируется не минимальным cut'ом, а «усреднённым» поведением.

**Множественные bottleneck'ы:** λ₃ − λ₂. Малая разность — несколько конкурирующих cut'ов с похожим conductance. Это признак, но не критерий: нужно смотреть локализацию собственных векторов.

**Локализация второго собственного вектора:**

    IPR₂ = Σ_i v_{2,i}⁴ / (Σ_i v_{2,i}²)²,     IPR₂^π = Σ_i w_i²,   w_i ∝ π_i · v_{2,i}²,

где π_i = d_i / vol(V). Интерпретация:

- делокализованный режим (модульность, макроскопический cut): IPR₂ ~ 1/n;
- локализованный режим (периферия): IPR₂ ~ 1/k, где k — размер носителя;
- оценка носителя: PR₂ = 1 / IPR₂.

Важно: IPR информативен только при **невырожденном** λ₂. Если λ₂ вырождено (например, у K_n кратность n−1), собственный вектор произволен и IPR не определён однозначно.

**Профиль conductance:** φ_curve(ρ) — полный sweep-профиль вдоль f₂ в обе стороны. Используется вместо фиксированного ε для определения режимов post hoc.

---

## 8. Фазовая диаграмма (гипотетическая)

```text
              макроскопический cut
                      ↑
                      │
       MODULAR         │       MIXED
                      │
                      │
──────────────────────┼──────────────────→ log Δ
                      │
      PERIPHERAL      │       CORE-DOMINATED
                      │
                      ↓
                локальный cut
```

---

## 9. Гипотезы

**H1.** Для sparse heavy-tailed configuration model с P(D = k) ∼ k^(−γ) асимптотика λ₂(G) определяется не только γ, но конкуренцией трёх масштабов:

    low-degree structures  ↔  community cuts  ↔  hub/core structure.

**H2.** Для полного графа при наличии положительной доли низкостепенных вершин λ₂(G) → 0 может быть тривиальным следствием периферии/несвязности.

**H3.** После перехода к giant component или 2-core λ₂(G_2core) может оставаться отделённой от нуля в режиме, где классическая Chung–Lu–Vu теория не применима.

**H4.** Если дополнительно наложить planted modular structure, то появляется другой масштаб:

    λ₂ ≍ φ_module

с точностью до квадратов/линейных факторов, задаваемых Cheeger inequality.

---

## 10. Что будет новым результатом

**Теорема (целевая):**

    λ₂(G_n) = Θ(f(n, γ))

в определённом режиме.

**Отдельно:**

    λ₂(G_n^2core) = Θ(g(γ, d̄))

при других условиях.

**Если точную асимптотику получить не удастся:**

    c·f(n, γ)  ≤  λ₂  ≤  C·f(n, γ)

— уже полноценный результат.

---

## 11. Математическая постановка

Конфигурационная модель с детерминированной последовательностью d₁, …, d_n:

    Σ_i d_i ≡ 0 (mod 2),

    N_k(n) = #{ i : d_i = k } ∼ c·n·k^(−γ).

Варьируем:

    d_min = 1, 2, 3, …

    d_max ∼ n^(1/(γ−1))

либо controlled cutoff.

**Методы:** message-passing подходы для произвольных degree distributions; spectral methods для конфигурационных моделей.

---

## 12. Экспериментальная часть

Для каждого (n, γ, d_min):

    n = 10³, 10⁴, 10⁵, …

генерируем ~100 реализаций.

Для каждой считаем:

    λ₂(G),   λ₂(G_giant),   λ₂(G_2core),

    φ(G),    φ(G_giant),    φ(G_2core),

    d_max,   d̄,   N₁,   N₂.

Scaling:

    log λ₂   vs   log n.

Если λ₂ ∼ n^(−α), то

    log λ₂ = −α·log n + C,

и α оценивается статистически.

**Четыре ансамбля:**

| Ансамбль | Что проверяем |
|---|---|
| Power-law d_min = 1 | периферия |
| Power-law d_min = 2 | цепочки / 2-core |
| Power-law d_min = 3 | экспандероподобный core |
| Power-law + planted communities | настоящая модульность |

---

## 13. Алгоритм вычисления φ_periph и φ_mod

**Подход 1: Спектральная аппроксимация через sweep cut**

1. Вычислить второй собственный вектор f₂ нормализованного лапласиана L = I − D^(−1/2)·A·D^(−1/2).
2. Отсортировать вершины по f₂(v): v₁, …, v_n.
3. Sweep cut: для каждого k = 1, …, n−1 рассмотреть S_k = {v₁, …, v_k} и вычислить φ(S_k).
4. Найти минимумы:
   - φ_periph^spec = min{ φ(S_k) : vol(S_k) ≤ ε·vol(V) }
   - φ_mod^spec = min{ φ(S_k) : ε·vol(V) ≤ vol(S_k) ≤ ½·vol(V) }
5. Соответствующие ρ_*.

**Сложность:** O(m·iter + n·log n).

**Подход 2: Local Cheeger с seed-множествами**

1. Выбрать seed-множества (низкостепенные вершины для периферии; хабы или случайные вершины для модульности).
2. Для каждого seed применить local clustering (Andersen–Chung–Lang, PageRank-based).
3. Собрать все найденные множества S, вычислить φ(S).
4. Отфильтровать по объёму.

**Сложность:** O(n·m·log n / φ²).

**Рекомендация:** начать со спектрального sweep cut. Если ρ₂ ≠ ρ_*, добавить local clustering для проверки.

---

## 14. Псевдокод sweep cut

```python
import numpy as np
from scipy.sparse.linalg import eigsh

def spectral_diagnostics(L, degrees, adj, eps=0.05):
    # L: нормализованный лапласиан
    # degrees: вектор степеней
    # adj: список смежности
    n = len(degrees)
    vol_V = degrees.sum()

    # Вычисляем lambda2 и f2
    vals, vecs = eigsh(L, k=2, which='SM')
    lambda2 = vals[1]
    f2 = vecs[:, 1]

    # Сортируем вершины по f2
    order = np.argsort(f2)

    vol_S = 0
    cut_size = 0
    in_S = np.zeros(n, dtype=bool)

    phi_periph = np.inf
    phi_mod = np.inf
    rho_periph = None
    rho_mod = None

    for k in range(n - 1):
        v = order[k]
        in_S[v] = True
        vol_S += degrees[v]

        for u in adj[v]:
            if in_S[u]:
                cut_size -= 1
            else:
                cut_size += 1

        if vol_S == 0:
            continue
        phi = cut_size / vol_S

        if vol_S <= eps * vol_V:
            if phi < phi_periph:
                phi_periph = phi
                rho_periph = vol_S / vol_V
        elif vol_S <= 0.5 * vol_V:
            if phi < phi_mod:
                phi_mod = phi
                rho_mod = vol_S / vol_V

    log_Delta = np.log(phi_mod) - np.log(phi_periph) if phi_periph > 0 else np.inf

    return {
        'lambda2': lambda2,
        'phi_periph': phi_periph,
        'phi_mod': phi_mod,
        'rho_periph': rho_periph,
        'rho_mod': rho_mod,
        'log_Delta': log_Delta,
        'rho_2': None  # вычислить отдельно
    }
```

Реализация в этом репозитории: `src/diagnostics.py` (см. раздел «Структура кода»).

---

## 15. Биологическая интерпретация

**Осторожность со словом «устойчивость».** λ₂ характеризует глобальную связность/expansion и через Cheeger связан с bottleneck'ами. Для random walk P = D^(−1)·A и μ₂ = 1 − λ₂. Малый λ₂ соответствует медленной релаксации случайного блуждания. Это можно интерпретировать как **структурную интегрированность сети**, но не автоматически как биологическую robustness.

**Нужно определить, какую именно устойчивость изучаем:**

- устойчивость к удалению узлов;
- устойчивость к удалению рёбер;
- сохранение giant component;
- устойчивость потоков/диффузии;
- устойчивость функциональных модулей.

**Масштаб функционального разделения через ρ_*:**

- ρ_* ≈ 0.5 — две большие половины (два полушария, два метаболических режима);
- ρ_* ≈ 0.1 — маленький модуль, отделённый от остального;
- ρ_* ≈ n^(−1) — один узел или маленькая группа.

**Формулировка для биологов:**

> В биологических сетях малая λ₂ сама по себе не идентифицирует функциональную модульность. Предлагаемый анализ дополнительно определяет масштаб минимального bottleneck и позволяет различать локальную периферию и макроскопическое разделение сети.

---

## 16. Связь с хабами

Тяжёлый хвост создаёт d_max ≫ d̄. Хабы сильно влияют на спектр adjacency matrix: необычно высокие степени могут создавать изолированные собственные значения с локализованными собственными векторами около хабов. Но это **не означает**, что хабы обязательно увеличивают λ₂.

**Важное различие:**

    λ_max(A)   и   λ₂(L)

чувствуют разные свойства сети. Хаб может резко изменить верхнюю часть спектра, практически не решив глобальный bottleneck.

---

## 17. Четырёхуровневая структура проекта

**Уровень I — математическая модель.**
Configuration model: (d₁, …, d_n), N_k(n) ∼ c·n·k^(−γ). Варьируем: γ, d_min, d_max, n.

**Уровень II — спектральная асимптотика.**
Исследуем: λ₂(G), λ₂(G_giant), λ₂(G_2core).

**Уровень III — механизм возникновения щели.**
Измеряем: φ_periph, φ_mod, Δ, ρ_*. Проверяем: λ₂ ↔ min(φ_periph, φ_mod).

**Уровень IV — биологическая интерпретация.**
Проверяем, соответствует ли наблюдаемая малая щель peripheral sparsity или macroscopic modular organization. И только после этого обсуждаем implications для robustness / functional integration.

---

## 18. Главный исследовательский вопрос

> Как низкостепенная периферия и тяжёлый хвост распределения степеней совместно определяют асимптотику второй собственной величины нормализованного лапласиана в разреженных случайных графах, и как с помощью пары (Δ, ρ_*) различить уменьшение λ₂, вызванное периферийным bottleneck, от уменьшения, вызванного макроскопической модульностью?

Потенциальный результат — не просто формула, а **теория режимов возникновения малой спектральной щели**.

---

## 19. Открытые вопросы и риски

- **Вычислимость φ:** точное вычисление NP-трудно; нужна явная аппроксимация (sweep cut или local clustering).
- **Устойчивость Δ:** использовать log Δ.
- **Связь с собственным вектором:** сравнить ρ₂ и ρ_*.
- **Множественные bottleneck'ы:** добавить λ₃ − λ₂.
- **2-core vs giant component:** уточнить, когда совпадают, когда различаются.
- **Размазанная периферия:** если low-degree узлы не компактны, φ_periph может быть неинформативен; проверить численно.
- **Аналитическая часть:** может оказаться неподъёмной; план B — численные гипотезы + частичные оценки.
- **Биологическая интерпретация:** требует осторожности; не называть λ₂ мерой «устойчивости» без уточнения.

---

## 20. Что дальше (исторический план — НЕ актуален)

> ⚠️ **Проект закрыт** (см. статус в начале и `LITERATURE.md`): результат
> воспроизводит известное (Samukhin–Dorogovtsev–Mendes, PRE 2008), препринт
> **не готовится**. Ниже — исходный план с отметками о статусе.

1. ✅ Генерация графов для четырёх ансамблей — сделано (`src/generators.py`).
2. ✅ Sweep cut для φ_periph, φ_mod, Δ, ρ_* — сделано (`src/diagnostics.py`).
3. ✅ Численные эксперименты n = 10³, 10⁴, 10⁵ — сделано (`results_merged`, `results_large`).
4. ✅ Фазовая диаграмма (log Δ, ρ_*) — сделано (`analyze_results.py`).
5. ✅ Гипотеза λ₂ ↔ min(φ_periph, φ_mod) — проверено; связь через Cheeger соблюдена.
6. ✅ Эмпирические гипотезы о порядке λ₂ — сделано; итог: **воспроизведение 2008**.
7. ⛔ Аналитический вывод — не требуется (результат известен, 2003–2008).
8. ⛔ Препринт (arXiv) — **не выпускается** (содержательной новизны нет).

---

## 21. Структура кода

```
graf/
├── README.md                  # этот документ (Unicode-обозначения)
├── paper.tex                  # вёрстка статьи (LaTeX, T2A/babel russian)
├── paper.pdf                  # скомпилированный препринт (RU)
├── paper_en.tex               # английская версия (для конференции)
├── paper_en.pdf               # скомпилированный препринт (EN)
├── requirements.txt
├── run_experiments.py         # оркестрация: 4 ансамбля × (n, γ, d_min), CSV-вывод
├── analyze_results.py         # α(n), α(φ), фазовая диаграмма, режимы, фигуры
├── export_graph.py            # экспорт сохранённых .npz -> edgelist/CSV/GraphML/GML/JSON
├── src/
│   ├── generators.py          # degree sequences + configuration model, 4 ансамбля
│   └── diagnostics.py         # λ₂ (full/giant/2-core), sweep cut, φ, Δ, ρ, IPR, φ_curve
└── tests/
    ├── test_spectral.py       # регрессионные тесты (sweep cut vs brute force, Cheeger, μ)
    └── test_known_spectra.py  # λ₂ на графах с известным спектром + IPR + φ_curve
```

Запуск:

```bash
pip install -r requirements.txt
python run_experiments.py --quick          # быстрый smoke test
python tests/test_spectral.py              # регрессионные тесты
python tests/test_known_spectra.py         # тесты на известных спектрах

# полный дизайн: малые n, несколько μ, фиксированный d_max, параллельно
python run_experiments.py --n 100 200 500 1000 3000 10000 \
    --gammas 2.1 2.5 3.0 --dmin 1 2 3 --mu 0.05 0.1 0.2 0.3 \
    --dmax 50 --reps 50 --jobs 8 --save-curves --out results/

# анализ: α(n), α(φ), фазовая диаграмма, режимы, фигуры
python analyze_results.py --csv results/diagnostics.csv --out analysis/

# сохранённые графы (.npz) -> читаемые форматы
python export_graph.py results/graphs/pl_dmin2_g2.5_n1000_d2_r0.npz \
    --format edgelist --out graph.txt
python export_graph.py results/graphs --format graphml --glob "planted_*" --out export/

pdflatex paper.tex && pdflatex paper.tex   # сборка PDF
```

---

## 22. Библиография

### Спектры случайных графов и нормализованный лапласиан

1. F. Chung. *Spectral Graph Theory*. CBMS Regional Conference Series in Mathematics, vol. 92, AMS, 1997.
2. F. Chung, L. Lu, V. Vu. **Spectra of random graphs with given expected degrees.** *Proceedings of the National Academy of Sciences* 100(11):6313–6318, 2003.
3. F. Chung, L. Lu. *Complex Graphs and Networks*. CBMS Regional Conference Series in Mathematics, vol. 107, AMS, 2006.
4. F. Chung, M. Radcliffe. **On the spectra of general random graphs.** *Electronic Journal of Combinatorics* 18(1), P215, 2011.
5. A. Coja-Oghlan, S. Lanka. **Finding planted partitions in random graphs with general degree distributions.** *SIAM Journal on Discrete Mathematics* 23(4):1682–1714, 2009.
6. A. Coja-Oghlan, S. Lanka. **The spectral gap of random graphs with given expected degrees.** *Electronic Journal of Combinatorics* 16(1), R138, 2009.

### Чигер, conductance и спектральная кластеризация

7. J. Cheeger. **A lower bound for the smallest eigenvalue of the Laplacian.** In *Problems in Analysis*, Princeton University Press, 195–199, 1970.
8. R. Kannan, S. Vempala, A. Vetta. **On clusterings: good, bad and spectral.** *Journal of the ACM* 51(3):497–515, 2004.
9. S. Arora, S. Rao, U. Vazirani. **Expander flows, geometric embeddings and graph partitioning.** *Journal of the ACM* 56(2), 2009.
10. R. Andersen, F. Chung, K. Lang. **Local graph partitioning using PageRank vectors.** *FOCS 2006*, 475–486.
11. J. R. Lee, S. Oveis Gharan, L. Trevisan. **Multi-way spectral partitioning and higher-order Cheeger inequalities.** *STOC 2012*, 1117–1130.

### Конфигурационная модель и структура разреженных графов

12. E. A. Bender, E. R. Canfield. **The asymptotic number of labeled graphs with given degree sequences.** *Journal of Combinatorial Theory, Series A* 24(3):296–307, 1978.
13. M. Molloy, B. Reed. **A critical point for random graphs with a given degree sequence.** *Random Structures & Algorithms* 6(2–3):161–180, 1995.
14. M. Molloy, B. Reed. **The size of the giant component of a random graph with a given degree sequence.** *Combinatorics, Probability and Computing* 7(3):295–305, 1998.
15. B. Pittel, J. Spencer, N. Wormald. **Sudden emergence of a giant k-core in a random graph.** *Journal of Combinatorial Theory, Series B* 67(1):111–151, 1996.
16. B. Bollobás, O. Riordan. **The diameter of a scale-free random graph.** *Combinatorica* 24(1):5–34, 2004.
17. S. Janson, T. Łuczak, A. Ruciński. *Random Graphs*. Wiley-Interscience, 2000.

### Степенное распределение, хабы и спектр сетей

18. A.-L. Barabási, R. Albert. **Emergence of scaling in random networks.** *Science* 286(5439):509–512, 1999.
19. M. Faloutsos, P. Faloutsos, C. Faloutsos. **On power-law relationships of the Internet topology.** *SIGCOMM 1999*, 251–262.
20. S. N. Dorogovtsev, A. V. Goltsev, J. F. F. Mendes, A. N. Samukhin. **Spectra of complex networks.** *Physical Review E* 68, 046109, 2003.
21. B. Söderberg. **General formalism for inhomogeneous random graphs.** *Physical Review E* 66, 066121, 2002.
22. C. Bordenave, M. Lelarge, L. Massoulié. **Non-backtracking spectrum of random graphs: community detection and non-regular Ramanujan graphs.** *Annals of Probability* 46(1):1–71, 2018.

### Модульность и сообщества

23. M. E. J. Newman, M. Girvan. **Finding and evaluating community structure in networks.** *Physical Review E* 69, 026113, 2004.
24. M. E. J. Newman. **Modularity and community structure in networks.** *PNAS* 103(23):8577–8582, 2006.
25. A. Decelle, F. Krzakala, C. Moore, L. Zdeborová. **Asymptotic analysis of the stochastic block model for modular networks and its algorithmic applications.** *Physical Review E* 84, 066106, 2011.
26. E. Abbe. **Community detection and stochastic block models: recent developments.** *Journal of Machine Learning Research* 18(177):1–86, 2018.

### Устойчивость биологических сетей

27. R. Albert, H. Jeong, A.-L. Barabási. **Error and attack tolerance of complex networks.** *Nature* 406:378–382, 2000.
28. R. Albert, A.-L. Barabási. **Statistical mechanics of complex networks.** *Reviews of Modern Physics* 74(1):47–97, 2002.
29. M. E. J. Newman. *Networks: An Introduction*. Oxford University Press, 2010.
30. A. Barrat, M. Barthélemy, A. Vespignani. *Dynamical Processes on Complex Networks*. Cambridge University Press, 2008.
