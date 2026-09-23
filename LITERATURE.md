# Литературная проверка новизны

Вопрос: публиковал ли кто-нибудь наблюдаемый скейлинг
`λ₂(normalized Laplacian, giant) ~ n^(−α)` при α ≈ 0.15…0.17
для power-law конфигурационной модели с **min degree 2** (d_max = 50).

Метод: arXiv API (полнотекстовые запросы) + Crossref. Semantic Scholar отдавал 429.
Дата: 2026-09-23.

---

## 🔴 ОБНОВЛЕНИЕ (после чтения ключевых работ): ядро УЖЕ ОПУБЛИКОВАНО

- **Samukhin, Dorogovtsev, Mendes**, *Laplacian spectra of complex networks and
  random walks on them: **Are scale-free architectures really important?***,
  **PRE 77, 036115 (2008)**, arXiv:0706.1176. Устанавливает: нижний край спектра
  `λ_c` определяется **минимальной степенью q_m**; «high-degree part … is not that
  essential»; **λ_c > 0 при q_m > 2, λ_c = 0 при q_m ≤ 2**. ⇒ наша дихотомия
  d_min=2/3 и вывод «хвост вторичен» — **воспроизведение известного результата**.
- **Dorogovtsev, Goltsev, Mendes, Samukhin**, *Spectra of complex networks*,
  **PRE 68, 046109 (2003)**: точные уравнения для спектра локально treelike сетей;
  `ρ(λ) ∼ |λ|^{1−2γ}`; отдельный разбор роли малых степеней.
- **Hata & Nakao**, *Localization of Laplacian eigenvectors on random networks*,
  **Sci Rep 7, 1121 (2017)**: локализация на вершинах близкой степени
  (degree–eigenvalue correspondence), но при `⟨k⟩≫1` (плотный режим).

**Вывод:** «новый закон» невозможен — ядро в литературе 2003–2008. Новый вопрос
должен быть **за** этим периметром (finite-size/флуктуации у перехода q_m=2→3;
локализация в разреженном режиме). Старые заметки ниже сохранены для истории.

---

## Итог (устарело, см. обновление выше)

**Точного совпадения** (именно `n^(−α)`, power-law config, min degree 2) **не
найдено** — предварительно ново. **Но** сам факт `λ₂ → 0` в разреженном режиме
уже известен (см. ER-прецедент), поэтому новизна возможна только в **скорости α
и её зависимости от распределения степеней**, а не в самом явлении.

## Что уже известно (база) — ВАЖНО для новизны

- **ISCCSP 2012**, DOI 10.1109/isccsp.2012.6217860: normalized Laplacian giant
  **разреженного ER** `G(n, d/n)` (d>1) имеет **нулевую щель** при n→∞.
  ⇒ качественный спад — **не** power-law-специфичен, а свойство разреженных
  графов с min degree ≤ 2 (листья + degree-2 цепочки).
- **Coja-Oghlan–Lanka**: у 2-core щель отделена от 0.
- **Hermon–Li–Yao–Zhang**: при min degree ≥ 3 relaxation time Θ(1).
  ⇒ переход d2→d3 = «выключение» degree-2-структуры.

## Что могло бы быть новым

Если при **равной средней степени** power-law спадает быстрее ER и это не
объясняется долей degree-2 — это и есть вклад тяжёлого хвоста. Проверяется
ER-контролем и контролями с ограниченным носителем (см. ниже).

## Ближайшие работы

| Работа | Что доказывает | Отношение к нам |
|---|---|---|
| **Hermon, Li, Yao, Zhang**, arXiv:2105.11585, *Ann. Prob.* 50(5) 2022 | config model со степенями в **[3, d̄]**, spectral gap / relaxation time | **ключевое**: требуют min degree ≥3 — ровно наш «const»-режим; d_min=2 вне условий |
| **Chung–Lu–Vu** (semicircle) | концентрация при **d̄_min ≫ ln²n** | наш режим d̄_min=O(1) — вне |
| **Wang, Li, Huang**, arXiv:2412.02087 (2024) | semicircle для норм. лапласиана config model, «mild assumptions» | про bulk-спектр, не про gap; heavy-tail+deg2, похоже, не покрыт |
| **Coste**, arXiv:1708.00530 (2017) | spectral gap разреженных random **digraph** по degree sequence (trace method) | методология |
| **Cai, Caputo, Perarnau, Quattropani**, arXiv:2104.08389 | directed CM, heavy-tailed in-degrees, mixing, power-law stationary | контекст heavy-tail |
| **Fernley**, arXiv:2211.13537 | voter model на scale-free, mixing, τ∈(2,3] | контекст «ultrasmall world» |
| **Lu & Peng**, arXiv:1204.6207 | Laplacian concentration требует δ ≫ ln n | почему sparse тяжёл |
| **Ducatez & Rivier**, EJP 2025, 10.1214/25-ejp1303 | gap (комбинаторного) лапласиана ER и вложенные деревья | методология для ER; не противоречит «нулевой щели» норм. лапласиана |
| **Avrachenkov et al.**, arXiv:1910.08869 | smallest nonzero eigenvalue → 0, power-law плотность у 0, spectral dimension | похожая феноменология |

## Что НЕ покрыто поиском (риски)

- Citation-граф: кто цитирует Coja-Oghlan–Lanka, Hermon et al., Bordenave — надёжнее всего, но arXiv API этого не даёт.
- Google Scholar / Web of Science / MathSciNet.
- Запросы про «spectral dimension» scale-free могли пропустить уже известный показатель.

## Наш ER-контроль (`run_er_control.py`, `results_er/`, 2026-09-23)

ER `G(n, c/n)`, c ∈ {2,3,4}, n = 100…10⁵, reps = 20. Сопоставление по средней
степени 2-core (ER-core плотнее всего графа):

| модель | core d̄ | α(λ₂_core) |
|---|---|---|
| ER c=2 | 2.68 | 0.138 |
| ER c=3 | 3.44 | 0.089 |
| ER c=4 | 4.26 | 0.072 |
| pl_dmin2 γ=3.0 | 3.08 | 0.171 |
| pl_dmin2 γ=2.5 | 3.87 | 0.158 |
| pl_dmin2 γ=2.1 | 5.06 | 0.150 |

**Вывод:** качественный спад генерический (ER тоже), но при равной d̄ у
power-law α **~2× больше**. Конфаунд: доля degree-2 у pl выше. Нужны контроли
с ограниченным носителем, чтобы отделить «хвост» от «доли degree-2».

## Что нужно для надёжной проверки новизны

1. **M_b/M_a контроли** (ограниченный носитель) — ключевое (см. CHANGELOG).
2. Targeted-поиск «кто цитирует» по Hermon–Li–Yao–Zhang / ISCCSP 2012.
3. Google Scholar / MathSciNet (не покрыто).
4. Тест `d_max`: варьировать длину хвоста при фиксированном n.
